#![allow(non_snake_case)]
use super::{ec_add_unequal, ec_select, ec_select_from_bits, EcPoint, EccChip};
use crate::halo2_proofs::{arithmetic::CurveAffine, group::Curve};
use halo2_base::{
    gates::{GateInstructions, RangeInstructions},
    halo2_proofs::{curves::pasta::pallas::Base, ff::PrimeField},
    utils::{fe_to_biguint, CurveAffineExt, ScalarField},
    AssignedValue, Context,
    QuantumCell::Existing,
};
use halo2_ecc::fields::FieldChip;
use itertools::Itertools;
use num_bigint::BigUint;
use std::{cmp::min, marker::PhantomData};

// this only works for curves GA with base field of prime order
#[derive(Clone, Debug)]
pub struct FixedEcPoint<BF: ScalarField> {
    pub x: BF,
    pub y: BF,
}

impl<BF: ScalarField> FixedEcPoint<BF> {
    pub fn construct(x: BF, y: BF) -> Self {
        Self { x, y }
    }

    pub fn from_curve<C: CurveAffineExt<Base = BF>>(point: C) -> Self {
        let (x, y) = point.into_coordinates();
        Self { x, y }
    }

    pub fn assign<'v, EC, Gate>(self, gate: &Gate, ctx: &mut Context<'_, BF>) -> EcPoint<BF>
    where
        Gate: GateInstructions<BF>,
    {
        let x = gate.load_constant(ctx, self.x);
        let y = gate.load_constant(ctx, self.y);
        EcPoint::construct(x, y)
    }

    pub fn assign_without_caching<'v>(self, ctx: &mut Context<'_, BF>) -> EcPoint<BF> {
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

pub fn scalar_multiply<'v, C, SFC>(
    ctx: &mut Context<'_, C::Base>,
    point: &C,
    scalar: &SFC::FieldPoint,
    max_bits: usize,
    window_bits: usize,
) -> EcPoint<C::Base>
where
    C: CurveAffineExt,
    C::Base: ScalarField,
    SFC: FieldChip<C::Base, FieldType = C::Scalar>,
{
    unimplemented!();
}

// basically just adding up individual fixed_base::scalar_multiply except that
// we do all batched normalization of cached points at once to further save
// inversion time during witness generation we also use the random accumulator
// for some extra efficiency (which also works in scalar multiply case but that
// is TODO)
pub fn msm<'v, C, SFC, Gate>(
    chip: &EccChip<C::Base, Gate>,
    ctx: &mut Context<'_, C::Base>,
    points: &[C],
    scalars: &[SFC::FieldPoint],
    max_scalar_bits_per_cell: usize,
    window_bits: usize,
) -> EcPoint<C::Base>
where
    C: CurveAffineExt,
    C::Base: ScalarField,
    SFC: FieldChip<C::Base, FieldType = C::Scalar>,
    Gate: GateInstructions<C::Base>,
{
    unimplemented!();
}
