use crate::halo2_verifier::{
    loader::{LoadedScalar, Loader},
    util::{
        arithmetic::{CurveAffine, Domain, Field, Fraction, Rotation},
        Itertools,
    },
    Protocol,
};
use num_integer::Integer;
use num_traits::One;
use serde::{Deserialize, Serialize};
use std::{
    cmp::max,
    collections::{BTreeMap, BTreeSet},
    fmt::Debug,
    iter::{self, Sum},
    ops::{Add, Mul, Neg, Sub},
};

impl<C> Protocol<C>
where
    C: CurveAffine,
{
    pub fn loaded<L: Loader<C>>(&self, loader: &L) -> Protocol<C, L> {
        let preprocessed = self
            .preprocessed
            .iter()
            .map(|preprocessed| loader.ec_point_load_const(preprocessed))
            .collect();
        let transcript_initial_state = self
            .transcript_initial_state
            .as_ref()
            .map(|transcript_initial_state| loader.load_const(transcript_initial_state));
        Protocol {
            domain: self.domain.clone(),
            preprocessed,
            num_instance: self.num_instance.clone(),
            num_witness: self.num_witness.clone(),
            num_challenge: self.num_challenge.clone(),
            evaluations: self.evaluations.clone(),
            queries: self.queries.clone(),
            quotient: self.quotient.clone(),
            transcript_initial_state,
            instance_committing_key: self.instance_committing_key.clone(),
            linearization: self.linearization,
            accumulator_indices: self.accumulator_indices.clone(),
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum CommonPolynomial {
    Identity,
    Lagrange(i32),
}

#[derive(Clone, Debug)]
pub struct CommonPolynomialEvaluation<C, L>
where
    C: CurveAffine,
    L: Loader<C>,
{
    zn: L::LoadedScalar,
    zn_minus_one: L::LoadedScalar,
    zn_minus_one_inv: Fraction<L::LoadedScalar>,
    identity: L::LoadedScalar,
    lagrange: BTreeMap<i32, Fraction<L::LoadedScalar>>,
}

impl<C, L> CommonPolynomialEvaluation<C, L>
where
    C: CurveAffine,
    L: Loader<C>,
{
    pub fn new(
        domain: &Domain<C::Scalar>,
        langranges: impl IntoIterator<Item = i32>,
        z: &L::LoadedScalar,
    ) -> Self {
        let loader = z.loader();

        let zn = z.pow_const(domain.n as u64);
        let langranges = langranges.into_iter().sorted().dedup().collect_vec();

        let one = loader.load_one();
        let zn_minus_one = zn.clone() - &one;
        let zn_minus_one_inv = Fraction::one_over(zn_minus_one.clone());

        let n_inv = loader.load_const(&domain.n_inv);
        let numer = zn_minus_one.clone() * &n_inv;
        let omegas = langranges
            .iter()
            .map(|&i| loader.load_const(&domain.rotate_scalar(C::Scalar::ONE, Rotation(i))))
            .collect_vec();
        let lagrange_evals = omegas
            .iter()
            .map(|omega| Fraction::new(numer.clone() * omega, z.clone() - omega))
            .collect_vec();

        Self {
            zn,
            zn_minus_one,
            zn_minus_one_inv,
            identity: z.clone(),
            lagrange: langranges.into_iter().zip(lagrange_evals).collect(),
        }
    }

    pub fn zn(&self) -> &L::LoadedScalar {
        &self.zn
    }

    pub fn zn_minus_one(&self) -> &L::LoadedScalar {
        &self.zn_minus_one
    }

    pub fn zn_minus_one_inv(&self) -> &L::LoadedScalar {
        self.zn_minus_one_inv.evaluated()
    }

    pub fn get(&self, poly: CommonPolynomial) -> &L::LoadedScalar {
        match poly {
            CommonPolynomial::Identity => &self.identity,
            CommonPolynomial::Lagrange(i) => self.lagrange.get(&i).unwrap().evaluated(),
        }
    }

    pub fn denoms(&mut self) -> impl IntoIterator<Item = &'_ mut L::LoadedScalar> {
        self.lagrange
            .iter_mut()
            .map(|(_, value)| value.denom_mut())
            .chain(iter::once(self.zn_minus_one_inv.denom_mut()))
            .flatten()
    }

    pub fn evaluate(&mut self) {
        self.lagrange
            .iter_mut()
            .map(|(_, value)| value)
            .chain(iter::once(&mut self.zn_minus_one_inv))
            .for_each(Fraction::evaluate)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QuotientPolynomial<F: Clone> {
    pub chunk_degree: usize,
    pub numerator: Expression<F>,
}

impl<F: Clone> QuotientPolynomial<F> {
    pub fn num_chunk(&self) -> usize {
        Integer::div_ceil(&(self.numerator.degree() - 1), &self.chunk_degree)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Deserialize, Serialize)]
pub struct Query {
    pub poly: usize,
    pub rotation: Rotation,
}

impl Query {
    pub fn new<R: Into<Rotation>>(poly: usize, rotation: R) -> Self {
        Self { poly, rotation: rotation.into() }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Expression<F> {
    Constant(F),
    CommonPolynomial(CommonPolynomial),
    Polynomial(Query),
    Challenge(usize),
    Negated(Box<Expression<F>>),
    Sum(Box<Expression<F>>, Box<Expression<F>>),
    Product(Box<Expression<F>>, Box<Expression<F>>),
    Scaled(Box<Expression<F>>, F),
    DistributePowers(Vec<Expression<F>>, Box<Expression<F>>),
}

impl<F: Clone> Expression<F> {
    pub fn evaluate<T: Clone>(
        &self,
        constant: &impl Fn(F) -> T,
        common_poly: &impl Fn(CommonPolynomial) -> T,
        poly: &impl Fn(Query) -> T,
        challenge: &impl Fn(usize) -> T,
        negated: &impl Fn(T) -> T,
        sum: &impl Fn(T, T) -> T,
        product: &impl Fn(T, T) -> T,
        scaled: &impl Fn(T, F) -> T,
    ) -> T {
        let evaluate = |expr: &Expression<F>| {
            expr.evaluate(constant, common_poly, poly, challenge, negated, sum, product, scaled)
        };
        match self {
            Expression::Constant(scalar) => constant(scalar.clone()),
            Expression::CommonPolynomial(poly) => common_poly(*poly),
            Expression::Polynomial(query) => poly(*query),
            Expression::Challenge(index) => challenge(*index),
            Expression::Negated(a) => {
                let a = evaluate(a);
                negated(a)
            }
            Expression::Sum(a, b) => {
                let a = evaluate(a);
                let b = evaluate(b);
                sum(a, b)
            }
            Expression::Product(a, b) => {
                let a = evaluate(a);
                let b = evaluate(b);
                product(a, b)
            }
            Expression::Scaled(a, scalar) => {
                let a = evaluate(a);
                scaled(a, scalar.clone())
            }
            Expression::DistributePowers(exprs, scalar) => {
                assert!(!exprs.is_empty());
                if exprs.len() == 1 {
                    return evaluate(exprs.first().unwrap());
                }
                let mut exprs = exprs.iter();
                let first = evaluate(exprs.next().unwrap());
                let scalar = evaluate(scalar);
                exprs.fold(first, |acc, expr| sum(product(acc, scalar.clone()), evaluate(expr)))
            }
        }
    }

    pub fn degree(&self) -> usize {
        match self {
            Expression::Constant(_) => 0,
            Expression::CommonPolynomial(_) => 1,
            Expression::Polynomial { .. } => 1,
            Expression::Challenge { .. } => 0,
            Expression::Negated(a) => a.degree(),
            Expression::Sum(a, b) => max(a.degree(), b.degree()),
            Expression::Product(a, b) => a.degree() + b.degree(),
            Expression::Scaled(a, _) => a.degree(),
            Expression::DistributePowers(a, b) => {
                a.iter().chain(Some(b.as_ref())).map(Self::degree).max().unwrap_or_default()
            }
        }
    }

    pub fn used_langrange(&self) -> BTreeSet<i32> {
        self.evaluate(
            &|_| None,
            &|poly| match poly {
                CommonPolynomial::Lagrange(i) => Some(BTreeSet::from_iter([i])),
                _ => None,
            },
            &|_| None,
            &|_| None,
            &|a| a,
            &merge_left_right,
            &merge_left_right,
            &|a, _| a,
        )
        .unwrap_or_default()
    }

    pub fn used_query(&self) -> BTreeSet<Query> {
        self.evaluate(
            &|_| None,
            &|_| None,
            &|query| Some(BTreeSet::from_iter([query])),
            &|_| None,
            &|a| a,
            &merge_left_right,
            &merge_left_right,
            &|a, _| a,
        )
        .unwrap_or_default()
    }
}

impl<F: Clone> From<Query> for Expression<F> {
    fn from(query: Query) -> Self {
        Self::Polynomial(query)
    }
}

impl<F: Clone> From<CommonPolynomial> for Expression<F> {
    fn from(common_poly: CommonPolynomial) -> Self {
        Self::CommonPolynomial(common_poly)
    }
}

macro_rules! impl_expression_ops {
    ($trait:ident, $op:ident, $variant:ident, $rhs:ty, $rhs_expr:expr) => {
        impl<F: Clone> $trait<$rhs> for Expression<F> {
            type Output = Expression<F>;
            fn $op(self, rhs: $rhs) -> Self::Output {
                Expression::$variant((self).into(), $rhs_expr(rhs).into())
            }
        }
        impl<F: Clone> $trait<$rhs> for &Expression<F> {
            type Output = Expression<F>;
            fn $op(self, rhs: $rhs) -> Self::Output {
                Expression::$variant((self.clone()).into(), $rhs_expr(rhs).into())
            }
        }
        impl<F: Clone> $trait<&$rhs> for Expression<F> {
            type Output = Expression<F>;
            fn $op(self, rhs: &$rhs) -> Self::Output {
                Expression::$variant((self).into(), $rhs_expr(rhs.clone()).into())
            }
        }
        impl<F: Clone> $trait<&$rhs> for &Expression<F> {
            type Output = Expression<F>;
            fn $op(self, rhs: &$rhs) -> Self::Output {
                Expression::$variant((self.clone()).into(), $rhs_expr(rhs.clone()).into())
            }
        }
    };
}

impl_expression_ops!(Mul, mul, Product, Expression<F>, std::convert::identity);
impl_expression_ops!(Mul, mul, Scaled, F, std::convert::identity);
impl_expression_ops!(Add, add, Sum, Expression<F>, std::convert::identity);
impl_expression_ops!(Sub, sub, Sum, Expression<F>, Neg::neg);

impl<F: Clone> Neg for Expression<F> {
    type Output = Expression<F>;
    fn neg(self) -> Self::Output {
        Expression::Negated(Box::new(self))
    }
}

impl<F: Clone> Neg for &Expression<F> {
    type Output = Expression<F>;
    fn neg(self) -> Self::Output {
        Expression::Negated(Box::new(self.clone()))
    }
}

impl<F: Clone + Default> Sum for Expression<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|acc, item| acc + item).unwrap_or_else(|| Expression::Constant(F::default()))
    }
}

impl<F: Field> One for Expression<F> {
    fn one() -> Self {
        Expression::Constant(F::ONE)
    }
}

fn merge_left_right<T: Ord>(a: Option<BTreeSet<T>>, b: Option<BTreeSet<T>>) -> Option<BTreeSet<T>> {
    match (a, b) {
        (Some(a), None) | (None, Some(a)) => Some(a),
        (Some(mut a), Some(b)) => {
            a.extend(b);
            Some(a)
        }
        _ => None,
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum LinearizationStrategy {
    /// Older linearization strategy of GWC19, which has linearization
    /// polynomial that doesn't evaluate to 0, and requires prover to send extra
    /// evaluation of it to verifier.
    WithoutConstant,
    /// Current linearization strategy of GWC19, which has linearization
    /// polynomial that evaluate to 0 by subtracting product of vanishing and
    /// quotient polynomials.
    MinusVanishingTimesQuotient,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct InstanceCommittingKey<C> {
    pub bases: Vec<C>,
    pub constant: Option<C>,
}
