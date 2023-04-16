use halo2_curves::group::ff::{Field, PrimeField};

use crate::backend::{
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::{inner_product, powers, product, BooleanHypercube},
        expression::{CommonPolynomial, Expression, Query},
        transcript::FieldTranscriptWrite,
        BitIndex, Itertools,
    },
    Error,
};
use std::{collections::HashMap, fmt::Debug, iter, sync::Arc};

pub mod classic;

#[derive(Clone, Debug)]
pub struct VirtualPolynomial<'a, F> {
    expression: &'a Expression<F>,
    polys: Vec<Arc<MultilinearPolynomial<F>>>,
    challenges: &'a [F],
    ys: &'a [Vec<F>],
}

impl<'a, F: PrimeField> VirtualPolynomial<'a, F> {
    pub fn new(
        expression: &'a Expression<F>,
        polys: impl IntoIterator<Item = Arc<MultilinearPolynomial<F>>>,
        challenges: &'a [F],
        ys: &'a [Vec<F>],
    ) -> Self {
        Self {
            expression,
            polys: polys.into_iter().collect(),
            challenges,
            ys,
        }
    }
}

pub trait SumCheck<F: Field>: Clone + Debug {
    type ProverParam: Clone + Debug;
    type VerifierParam: Clone + Debug;
    type Proof: Clone + Debug;

    fn prove(
        pp: &Self::ProverParam,
        num_vars: usize,
        virtual_poly: VirtualPolynomial<F>,
        sum: F,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(Vec<Vec<F>>, Vec<F>, Vec<F>), Error>;

    fn verify(
        vp: &Self::VerifierParam,
        msgs: &[&[F]],
        num_vars: usize,
        degree: usize,
        sum: F,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(F, Vec<F>), Error>;
}

pub fn evaluate<F: PrimeField>(
    expression: &Expression<F>,
    num_vars: usize,
    evals: &HashMap<Query, F>,
    challenges: &[F],
    ys: &[&[F]],
    x: &[F],
) -> F {
    assert!(num_vars > 0 && expression.max_used_rotation_distance() <= num_vars);
    let lagranges = {
        // let bh = BooleanHypercube::new(num_vars).iter().collect_vec();
        expression
            .used_langrange()
            .into_iter()
            .map(|i| {
                let b = i.rem_euclid(1 << num_vars as i32) as usize;
                (i, lagrange_eval(x, b))
            })
            .collect::<HashMap<_, _>>()
    };
    let eq_xys = ys.iter().map(|y| eq_xy_eval(x, y)).collect_vec();
    let eq_xy_shifted_00s = ys
        .iter()
        .map(|y| eq_xy_eval_shift_ab(x, y, false, false))
        .collect_vec();
    let eq_xy_shifted_01s = ys
        .iter()
        .map(|y| eq_xy_eval_shift_ab(x, y, false, true))
        .collect_vec();
    let eq_xy_shifted_10s = ys
        .iter()
        .map(|y| eq_xy_eval_shift_ab(x, y, true, false))
        .collect_vec();
    let eq_xy_shifted_11s = ys
        .iter()
        .map(|y| eq_xy_eval_shift_ab(x, y, true, true))
        .collect_vec();
    let identity = identity_eval(x);
    expression.evaluate(
        &|scalar| scalar,
        &|poly| match poly {
            CommonPolynomial::Lagrange(i) => lagranges[&i],
            CommonPolynomial::EqXY(idx) => eq_xys[idx],
            CommonPolynomial::EqXYShifted00(idx) => eq_xy_shifted_00s[idx],
            CommonPolynomial::EqXYShifted01(idx) => eq_xy_shifted_01s[idx],
            CommonPolynomial::EqXYShifted10(idx) => eq_xy_shifted_10s[idx],
            CommonPolynomial::EqXYShifted11(idx) => eq_xy_shifted_11s[idx],
            CommonPolynomial::Identity(idx) => F::from((idx << num_vars) as u64) + identity,
        },
        &|query| evals[&query],
        &|idx| challenges[idx],
        &|scalar| -scalar,
        &|lhs, rhs| lhs + &rhs,
        &|lhs, rhs| lhs * &rhs,
        &|value, scalar| scalar * value,
    )
}

pub fn lagrange_eval<F: PrimeField>(x: &[F], b: usize) -> F {
    assert!(!x.is_empty());

    product(x.iter().enumerate().map(
        |(idx, x_i)| {
            if b.nth_bit(idx) {
                *x_i
            } else {
                F::one() - x_i
            }
        },
    ))
}

pub fn eq_xy_eval<F: PrimeField>(x: &[F], y: &[F]) -> F {
    assert!(!x.is_empty());
    assert_eq!(x.len(), y.len());

    product(
        x.iter()
            .zip(y)
            .map(|(x_i, y_i)| (*x_i * y_i).double() + F::one() - x_i - y_i),
    )
}

/// Compute ((1 - yn) * (1 - a) + yn * a ) eq ((b, y1, ..., yn), x)
pub fn eq_xy_eval_shift_ab<F: PrimeField>(x: &[F], y: &[F], a: bool, b: bool) -> F {
    assert!(!x.is_empty());
    let y_last = if a {
        y[y.len() - 1]
    } else {
        F::one() - y[y.len() - 1]
    };
    let y_0 = if b { F::one() } else { F::zero() };
    let y = iter::once(y_0)
        .chain(y[..y.len() - 1].to_vec())
        .collect_vec();
    y_last * eq_xy_eval(x, &y)
}

fn identity_eval<F: PrimeField>(x: &[F]) -> F {
    inner_product(x, &powers(F::from(2)).take(x.len()).collect_vec())
}

#[cfg(test)]
pub(super) mod test {
    use crate::backend::{
        piop::sum_check::{evaluate, SumCheck, VirtualPolynomial},
        poly::multilinear::{rotation_eval, MultilinearPolynomial},
        util::{expression::Expression, transcript::Keccak256Transcript},
    };
    use halo2_curves::bn256::Fr;
    use itertools::Itertools;
    use std::{ops::Range, sync::Arc};

    pub fn run_sum_check<S: SumCheck<Fr>>(
        num_vars: usize,
        expression_fn: impl Fn(usize) -> Expression<Fr>,
        param_fn: impl Fn(usize) -> (S::ProverParam, S::VerifierParam),
        assignment_fn: impl Fn(usize) -> (Vec<Arc<MultilinearPolynomial<Fr>>>, Vec<Fr>, Vec<Fr>),
        sum: Fr,
    ) -> (Vec<Fr>, Vec<Fr>) {
        let expression = expression_fn(num_vars);
        let degree = expression.degree();
        let (pp, vp) = param_fn(expression.degree());
        let (polys, challenges, y) = assignment_fn(num_vars);

        let ys = [y];
        let (proof, point, poly_evals) = {
            let virtual_poly =
                VirtualPolynomial::new(&expression, polys.iter().cloned(), &challenges, &ys);
            let mut transcript = Keccak256Transcript::<Vec<u8>>::default();
            S::prove(&pp, num_vars, virtual_poly, sum, &mut transcript).unwrap()
        };
        let accept = {
            let mut transcript = Keccak256Transcript::<Vec<u8>>::default();
            let proof_slices = proof.iter().map(|p| p.as_slice()).collect_vec();
            let (x_eval, x) =
                S::verify(&vp, &proof_slices, num_vars, degree, sum, &mut transcript).unwrap();
            let evals = expression
                .used_query()
                .into_iter()
                .map(|query| {
                    let evaluate_for_rotation =
                        polys[query.poly()].evaluate_for_rotation(&x, query.rotation());
                    let eval = rotation_eval(&x, query.rotation(), &evaluate_for_rotation);
                    (query, eval)
                })
                .collect();
            x_eval == evaluate(&expression, num_vars, &evals, &challenges, &[&ys[0]], &x)
        };
        assert!(accept);
        (point, poly_evals)
    }

    pub fn run_zero_check<S: SumCheck<Fr>>(
        num_vars_range: Range<usize>,
        expression_fn: impl Fn(usize) -> Expression<Fr>,
        param_fn: impl Fn(usize) -> (S::ProverParam, S::VerifierParam),
        assignment_fn: impl Fn(usize) -> (Vec<Arc<MultilinearPolynomial<Fr>>>, Vec<Fr>, Vec<Fr>),
    ) {
        for num_vars in num_vars_range {
            let expression = expression_fn(num_vars);
            let degree = expression.degree();
            run_sum_check::<S>(
                num_vars,
                |_| expression.clone(),
                |_| param_fn(degree),
                |_| assignment_fn(num_vars),
                Fr::zero(),
            );
        }
    }

    pub fn run_constraint_opening_check<S: SumCheck<Fr>>(
        num_vars_range: Range<usize>,
        expression_fn: impl Fn(usize) -> (Expression<Fr>, Expression<Fr>),
        param_fn: impl Fn(usize) -> (S::ProverParam, S::VerifierParam),
        assignment_fn: impl Fn(usize) -> (Vec<Arc<MultilinearPolynomial<Fr>>>, Vec<Fr>, Vec<Fr>),
    ) {
        for num_vars in num_vars_range {
            let (constraint_expression, opening_expression) = expression_fn(num_vars);
            let constraint_degree = constraint_expression.degree();
            let param = param_fn(constraint_degree);
            let (poly, challenges, y) = assignment_fn(num_vars);
            let eta = challenges[challenges.len() - 1];
            println!("Sumcheck for constraints");
            let (point, evals) = run_sum_check::<S>(
                num_vars,
                |_| constraint_expression.clone(),
                |_| param.clone(),
                |_| (poly.clone(), challenges.clone(), y.clone()),
                Fr::zero(),
            );
            println!("Sumcheck for openings");
            let sum = evals
                .iter()
                .skip(1) // skip instance_polys
                .fold(Fr::zero(), |acc, eval| acc * eta + eval);
            run_sum_check::<S>(
                num_vars,
                |_| opening_expression.clone(),
                |_| param.clone(),
                |_| (poly.clone(), challenges.clone(), point.clone()),
                sum,
            );
        }
    }

    macro_rules! tests {
        ($impl:ty) => {
            // #[test]
            // fn sum_check_lagrange() {
            //     use halo2_curves::bn256::Fr;
            //     use rand::rngs::OsRng;
            //     use std::sync::Arc;
            //     use $crate::backend::{
            //         piop::sum_check::test::run_zero_check,
            //         poly::multilinear::MultilinearPolynomial,
            //         util::{
            //             arithmetic::{BooleanHypercube, Field},
            //             expression::{CommonPolynomial, Expression, Query, Rotation},
            //             test::rand_vec,
            //             Itertools,
            //         },
            //     };

            //     run_zero_check::<$impl>(
            //         2..4,
            //         |num_vars| {
            //             let polys = (0..1 << num_vars)
            //                 .map(|idx| {
            //                     Expression::<Fr>::Polynomial(Query::new(idx, Rotation::cur()))
            //                 })
            //                 .collect_vec();
            //             let gates = polys
            //                 .iter()
            //                 .enumerate()
            //                 .map(|(i, poly)| {
            //                     Expression::CommonPolynomial(CommonPolynomial::Lagrange(i as i32))
            //                         - poly
            //                 })
            //                 .collect_vec();
            //             let alpha = Expression::Challenge(0);
            //             let eq = Expression::eq_xy(0);
            //             Expression::distribute_powers(&gates, &alpha) * eq
            //         },
            //         |_| ((), ()),
            //         |num_vars| {
            //             let polys = BooleanHypercube::new(num_vars)
            //                 .iter()
            //                 .map(|idx| {
            //                     let mut polys =
            //                         MultilinearPolynomial::new(vec![Fr::zero(); 1 << num_vars]);
            //                     polys[idx] = Fr::one();
            //                     Arc::new(polys)
            //                 })
            //                 .collect_vec();
            //             let alpha = Fr::random(OsRng);
            //             (polys, vec![alpha], rand_vec(num_vars, OsRng))
            //         },
            //     );
            // }

            #[test]
            fn sum_check_rotation() {
                use halo2_curves::bn256::Fr;
                use rand::rngs::OsRng;
                use std::{iter, sync::Arc};
                use $crate::backend::{
                    piop::sum_check::test::run_zero_check,
                    poly::multilinear::MultilinearPolynomial,
                    util::{
                        arithmetic::{BooleanHypercube, Field},
                        expression::{Expression, Query, Rotation},
                        test::rand_vec,
                        Itertools,
                    },
                };

                run_zero_check::<$impl>(
                    2..16,
                    |num_vars| {
                        let polys = (-(num_vars as i32) + 1..num_vars as i32)
                            .rev()
                            .enumerate()
                            .map(|(idx, rotation)| {
                                Expression::<Fr>::Polynomial(Query::new(idx, rotation.into()))
                            })
                            .collect_vec();
                        let gates = polys
                            .windows(2)
                            .map(|polys| &polys[1] - &polys[0])
                            .collect_vec();
                        let alpha = Expression::Challenge(0);
                        let eq = Expression::eq_xy(0);
                        Expression::distribute_powers(&gates, &alpha) * eq
                    },
                    |_| ((), ()),
                    |num_vars| {
                        let bh = BooleanHypercube::new(num_vars);
                        let rotate = |f: &Vec<Fr>| {
                            (0..1 << num_vars)
                                .map(|idx| f[bh.rotate(idx, Rotation::next())])
                                .collect_vec()
                        };
                        let poly = rand_vec(1 << num_vars, OsRng);
                        let polys = iter::successors(Some(poly), |poly| Some(rotate(poly)))
                            .map(|p| Arc::new(MultilinearPolynomial::new(p)))
                            .take(2 * num_vars - 1)
                            .collect_vec();
                        let alpha = Fr::random(OsRng);
                        (polys, vec![alpha], rand_vec(num_vars, OsRng))
                    },
                );
            }

            #[test]
            fn sum_check_plonk() {
                use halo2_curves::bn256::Fr;
                use rand::rngs::OsRng;
                use $crate::backend::{
                    piop::{
                        sum_check::test::run_constraint_opening_check,
                        util::{plonk_expression, rand_plonk_assignment},
                    },
                    util::test::rand_vec,
                };

                run_constraint_opening_check::<$impl>(
                    2..16,
                    |_| plonk_expression(),
                    |_| ((), ()),
                    |num_vars| {
                        let (polys, challenges) = rand_plonk_assignment(num_vars, OsRng);
                        (polys, challenges, rand_vec(num_vars, OsRng))
                    },
                );
            }

            #[test]
            fn sum_check_plonk_with_lookup() {
                use halo2_curves::bn256::Fr;
                use rand::rngs::OsRng;
                use $crate::backend::{
                    piop::{
                        sum_check::test::run_constraint_opening_check,
                        util::{plonk_with_lookup_expression, rand_plonk_with_lookup_assignment},
                    },
                    util::test::rand_vec,
                };

                run_constraint_opening_check::<$impl>(
                    2..16,
                    |_| plonk_with_lookup_expression(),
                    |_| ((), ()),
                    |num_vars| {
                        let (polys, challenges) =
                            rand_plonk_with_lookup_assignment(num_vars, OsRng);
                        (polys, challenges, rand_vec(num_vars, OsRng))
                    },
                );
            }
        };
    }

    pub(super) use tests;
}
