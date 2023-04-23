#![allow(non_snake_case)]
use super::{ec_add_unequal, ec_select, ec_select_from_bits, EcPoint, EccChip};
use crate::halo2_proofs::{arithmetic::CurveAffine, group::Curve};
use halo2_base::{
    gates::{GateInstructions, RangeInstructions},
    utils::{fe_to_biguint, CurveAffineExt, PrimeField},
    AssignedValue, Context,
    QuantumCell::Existing,
};
use itertools::Itertools;
use num_bigint::BigUint;
use std::{cmp::min, marker::PhantomData};

// this only works for curves GA with base field of prime order
#[derive(Clone, Debug)]
pub struct FixedEcPoint<F: PrimeField, C: CurveAffine> {
    pub x: AssignedValue<F>,
    pub y: AssignedValue<F>,
    _marker: PhantomData<C>,
}

impl<F: PrimeField, C: CurveAffineExt> FixedEcPoint<F, C>
where
    C::Base: PrimeField,
{
    pub fn construct(x: AssignedValue<F>, y: AssignedValue<F>) -> Self {
        Self {
            x,
            y,
            _marker: PhantomData,
        }
    }

    pub fn from_curve(point: C) -> Self {
        unimplemented!();
    }

    pub fn assign<'v, FC>(self, ctx: &mut Context<'_, F>) -> EcPoint<F> {
        unimplemented!();
    }

    pub fn assign_without_caching<'v, FC>(self, ctx: &mut Context<'_, F>) -> EcPoint<F> {
        unimplemented!();
    }
}

// computes `[scalar] * P` on y^2 = x^3 + b where `P` is fixed (constant)
// - `scalar` is represented as a reference array of `AssignedCell`s
// - `scalar = sum_i scalar_i * 2^{max_bits * i}`
// - an array of length > 1 is needed when `scalar` exceeds the modulus of
//   scalar field `F`
// assumes:
// - `scalar_i < 2^{max_bits} for all i` (constrained by num_to_bits)
// - `max_bits <= modulus::<F>.bits()`

pub fn scalar_multiply<'v, F, C>(
    ctx: &mut Context<'_, F>,
    point: &C,
    scalar: &[AssignedValue<F>],
    max_bits: usize,
    window_bits: usize,
) -> EcPoint<F>
where
    F: PrimeField,
    C: CurveAffineExt,
    C::Base: PrimeField,
{
    unimplemented!();
}

// basically just adding up individual fixed_base::scalar_multiply except that
// we do all batched normalization of cached points at once to further save
// inversion time during witness generation we also use the random accumulator
// for some extra efficiency (which also works in scalar multiply case but that
// is TODO)
pub fn msm<'v, F, C>(
    chip: &EccChip<F>,
    ctx: &mut Context<'_, F>,
    points: &[C],
    scalars: &[Vec<AssignedValue<F>>],
    max_scalar_bits_per_cell: usize,
    window_bits: usize,
) -> EcPoint<F>
where
    F: PrimeField,
    C: CurveAffineExt,
    C::Base: PrimeField,
{
    unimplemented!();
}
