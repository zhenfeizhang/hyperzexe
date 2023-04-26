use crate::hyperplonk_verifier::util::Itertools;
use num_bigint::BigUint;
use num_traits::One;
use serde::{Deserialize, Serialize};
use std::{
    cmp::Ordering,
    fmt::Debug,
    iter, mem,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::hyperplonk_verifier::halo2_curves;
pub use halo2_curves::{
    group::{
        ff::{BatchInvert, Field, PrimeField},
        prime::PrimeCurveAffine,
        Curve, Group, GroupEncoding,
    },
    pairing::MillerLoopResult,
    Coordinates, CurveAffine,
};

pub trait MultiMillerLoop: halo2_curves::pairing::MultiMillerLoop + Debug {}

impl<M: halo2_curves::pairing::MultiMillerLoop + Debug> MultiMillerLoop for M {}

pub trait FieldOps:
    Sized
    + Neg<Output = Self>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
    + for<'a> MulAssign<&'a Self>
{
    fn invert(&self) -> Option<Self>;
}

pub fn batch_invert_and_mul<F: PrimeField>(values: &mut [F], coeff: &F) {
    let products = values
        .iter()
        .filter(|value| !value.is_zero_vartime())
        .scan(F::ONE, |acc, value| {
            *acc *= value;
            Some(*acc)
        })
        .collect_vec();

    let mut all_product_inv = products.last().unwrap().invert().unwrap() * coeff;

    for (value, product) in values
        .iter_mut()
        .rev()
        .filter(|value| !value.is_zero_vartime())
        .zip(products.into_iter().rev().skip(1).chain(Some(F::ONE)))
    {
        let mut inv = all_product_inv * product;
        mem::swap(value, &mut inv);
        all_product_inv *= inv;
    }
}

pub fn batch_invert<F: PrimeField>(values: &mut [F]) {
    batch_invert_and_mul(values, &F::ONE)
}

pub fn root_of_unity<F: PrimeField>(k: usize) -> F {
    assert!(k <= F::S as usize);

    iter::successors(Some(F::ROOT_OF_UNITY), |acc| Some(acc.square()))
        .take(F::S as usize - k + 1)
        .last()
        .unwrap()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Rotation(pub i32);

impl Rotation {
    pub fn cur() -> Self {
        Rotation(0)
    }

    pub fn prev() -> Self {
        Rotation(-1)
    }

    pub fn next() -> Self {
        Rotation(1)
    }
}

impl From<i32> for Rotation {
    fn from(rotation: i32) -> Self {
        Self(rotation)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Domain<F: PrimeField> {
    pub k: usize,
    pub n: usize,
    pub n_inv: F,
    pub gen: F,
    pub gen_inv: F,
}

impl<F: PrimeField> Domain<F> {
    pub fn new(k: usize, gen: F) -> Self {
        let n = 1 << k;
        let n_inv = F::from(n as u64).invert().unwrap();
        let gen_inv = gen.invert().unwrap();

        Self { k, n, n_inv, gen, gen_inv }
    }

    pub fn rotate_scalar(&self, scalar: F, rotation: Rotation) -> F {
        match rotation.0.cmp(&0) {
            Ordering::Equal => scalar,
            Ordering::Greater => scalar * self.gen.pow_vartime([rotation.0 as u64]),
            Ordering::Less => scalar * self.gen_inv.pow_vartime([(-rotation.0) as u64]),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Fraction<T> {
    numer: Option<T>,
    denom: T,
    eval: Option<T>,
    inv: bool,
}

impl<T> Fraction<T> {
    pub fn new(numer: T, denom: T) -> Self {
        Self { numer: Some(numer), denom, eval: None, inv: false }
    }

    pub fn one_over(denom: T) -> Self {
        Self { numer: None, denom, eval: None, inv: false }
    }

    pub fn denom(&self) -> Option<&T> {
        if !self.inv {
            Some(&self.denom)
        } else {
            None
        }
    }

    pub fn denom_mut(&mut self) -> Option<&mut T> {
        if !self.inv {
            self.inv = true;
            Some(&mut self.denom)
        } else {
            None
        }
    }
}

impl<T: FieldOps + Clone> Fraction<T> {
    pub fn evaluate(&mut self) {
        assert!(self.inv);
        assert!(self.eval.is_none());

        self.eval = Some(
            self.numer
                .take()
                .map(|numer| numer * &self.denom)
                .unwrap_or_else(|| self.denom.clone()),
        );
    }

    pub fn evaluated(&self) -> &T {
        assert!(self.eval.is_some());

        self.eval.as_ref().unwrap()
    }
}

pub fn ilog2(value: usize) -> usize {
    (usize::BITS - value.leading_zeros() - 1) as usize
}

pub fn modulus<F: PrimeField>() -> BigUint {
    fe_to_big(-F::ONE) + 1usize
}

pub fn fe_from_big<F: PrimeField>(big: BigUint) -> F {
    let bytes = big.to_bytes_le();
    let mut repr = F::Repr::default();
    assert!(bytes.len() <= repr.as_ref().len());
    repr.as_mut()[..bytes.len()].clone_from_slice(bytes.as_slice());
    F::from_repr(repr).unwrap()
}

pub fn fe_to_big<F: PrimeField>(fe: F) -> BigUint {
    BigUint::from_bytes_le(fe.to_repr().as_ref())
}

pub fn fe_to_fe<F1: PrimeField, F2: PrimeField>(fe: F1) -> F2 {
    fe_from_big(fe_to_big(fe) % modulus::<F2>())
}

pub fn fe_from_limbs<F1: PrimeField, F2: PrimeField, const LIMBS: usize, const BITS: usize>(
    limbs: [F1; LIMBS],
) -> F2 {
    fe_from_big(
        limbs
            .iter()
            .map(|limb| BigUint::from_bytes_le(limb.to_repr().as_ref()))
            .zip((0usize..).step_by(BITS))
            .map(|(limb, shift)| limb << shift)
            .reduce(|acc, shifted| acc + shifted)
            .unwrap(),
    )
}

pub fn fe_to_limbs<F1: PrimeField, F2: PrimeField, const LIMBS: usize, const BITS: usize>(
    fe: F1,
) -> [F2; LIMBS] {
    let big = BigUint::from_bytes_le(fe.to_repr().as_ref());
    let mask = &((BigUint::one() << BITS) - 1usize);
    (0usize..)
        .step_by(BITS)
        .take(LIMBS)
        .map(|shift| fe_from_big((&big >> shift) & mask))
        .collect_vec()
        .try_into()
        .unwrap()
}

pub fn powers<F: Field>(scalar: F) -> impl Iterator<Item = F> {
    iter::successors(Some(F::ONE), move |power| Some(scalar * power))
}

pub fn inner_product<F: Field>(lhs: &[F], rhs: &[F]) -> F {
    lhs.iter()
        .zip_eq(rhs.iter())
        .map(|(lhs, rhs)| *lhs * rhs)
        .reduce(|acc, product| acc + product)
        .unwrap_or_default()
}
