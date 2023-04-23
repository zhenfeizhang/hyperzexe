#![allow(non_snake_case)]
use crate::{
    fields::{
        fp::{FpConfig, ScalarFieldChip},
        FieldChip, PrimeFieldChip, Selectable,
    },
    halo2_proofs::{
        arithmetic::CurveAffine,
        circuit::Value,
        group::{Curve, Group},
    },
};
use halo2_base::{
    gates::{flex_gate::FlexGateConfig, GateInstructions, RangeInstructions},
    halo2_proofs::group::prime::PrimeCurveAffine,
    utils::{modulus, CurveAffineExt, PrimeField},
    AssignedValue, Context,
    QuantumCell::Existing,
};
use itertools::Itertools;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::marker::PhantomData;

pub mod fixed_base;
// pub mod fixed_base_pippenger;
pub mod pippenger;

// EcPoint and EccChip, native implementation
#[derive(Debug)]
pub struct EcPoint<BF: PrimeField> {
    pub x: AssignedValue<BF>,
    pub y: AssignedValue<BF>,
}

impl<BF: PrimeField> Clone for EcPoint<BF> {
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.clone(),
        }
    }
}

impl<BF: PrimeField> EcPoint<BF> {
    pub fn construct(x: AssignedValue<BF>, y: AssignedValue<BF>) -> Self {
        Self { x, y }
    }

    pub fn x(&self) -> &AssignedValue<BF> {
        &self.x
    }

    pub fn y(&self) -> &AssignedValue<BF> {
        &self.y
    }
}

// Implements:
//  Given P = (x_1, y_1) and Q = (x_2, y_2), ecc points over the field F_p
//      assume x_1 != x_2
//  Find ec addition P + Q = (x_3, y_3)
// By solving:
//  lambda = (y_2-y_1)/(x_2-x_1) using constraint
//  lambda * (x_2 - x_1) = y_2 - y_1
//  x_3 = lambda^2 - x_1 - x_2 (mod p)
//  y_3 = lambda (x_1 - x_3) - y_1 mod p
//
/// For optimization reasons, we assume that if you are using this with
/// `is_strict = true`, then you have already called `chip.enforce_less_than_p`
/// on both `P.x` and `P.y`
pub fn ec_add_unequal<'v, BF: PrimeField>(
    ctx: &mut Context<'_, BF>,
    P: &EcPoint<BF>,
    Q: &EcPoint<BF>,
    is_strict: bool,
) -> EcPoint<BF> {
    unimplemented!();
}

// Implements:
//  Given P = (x_1, y_1) and Q = (x_2, y_2), ecc points over the field F_p
//  Find ecc subtraction P - Q = (x_3, y_3)
//  -Q = (x_2, -y_2)
//  lambda = -(y_2+y_1)/(x_2-x_1) using constraint
//  x_3 = lambda^2 - x_1 - x_2 (mod p)
//  y_3 = lambda (x_1 - x_3) - y_1 mod p
//  Assumes that P !=Q and Q != (P - Q)
//
/// For optimization reasons, we assume that if you are using this with
/// `is_strict = true`, then you have already called `chip.enforce_less_than_p`
/// on both `P.x` and `P.y`
pub fn ec_sub_unequal<'v, BF: PrimeField>(
    ctx: &mut Context<'_, BF>,
    P: &EcPoint<BF>,
    Q: &EcPoint<BF>,
    is_strict: bool,
) -> EcPoint<BF> {
    unimplemented!();
}

// Implements:
// computing 2P on elliptic curve E for P = (x, y)
// formula from https://crypto.stanford.edu/pbc/notes/elliptic/explicit.html
// assume y != 0 (otherwise 2P = O)

// lamb =  3x^2 / (2 y) % p
// x_3 = out[0] = lambda^2 - 2 x % p
// y_3 = out[1] = lambda (x - x_3) - y % p

// we precompute lambda and constrain (2y) * lambda = 3 x^2 (mod p)
// then we compute x_3 = lambda^2 - 2 x (mod p)
//                 y_3 = lambda (x - x_3) - y (mod p)
pub fn ec_double<'v, BF: PrimeField>(ctx: &mut Context<'_, BF>, P: &EcPoint<BF>) -> EcPoint<BF> {
    unimplemented!();
}

pub fn ec_select<'v, BF: PrimeField>(
    ctx: &mut Context<'_, BF>,
    P: &EcPoint<BF>,
    Q: &EcPoint<BF>,
    sel: &AssignedValue<BF>,
) -> EcPoint<BF> {
    unimplemented!();
}

// takes the dot product of points with sel, where each is intepreted as
// a _vector_
pub fn ec_select_by_indicator<'v, BF: PrimeField>(
    ctx: &mut Context<'_, BF>,
    points: &[EcPoint<BF>],
    coeffs: &[AssignedValue<BF>],
) -> EcPoint<BF> {
    unimplemented!();
}

// `sel` is little-endian binary
pub fn ec_select_from_bits<'v, BF: PrimeField>(
    ctx: &mut Context<'_, BF>,
    points: &[EcPoint<BF>],
    sel: &[AssignedValue<BF>],
) -> EcPoint<BF> {
    unimplemented!();
}

// computes [scalar] * P on y^2 = x^3 + b
//   * P has order given by the scalar field modulus
pub fn scalar_multiply<'v, BF: PrimeField, C>(
    ctx: &mut Context<'_, BF>,
    P: &EcPoint<BF>,
    scalar: &<ScalarFieldChip<C> as FieldChip<BF>>::FieldPoint,
    max_bits: usize,
    window_bits: usize,
) -> EcPoint<BF>
where
    C: CurveAffineExt<Base = BF>,
{
    unimplemented!();
}

pub fn is_on_curve<'v, BF, C>(ctx: &mut Context<'_, BF>, P: &EcPoint<BF>)
where
    BF: PrimeField,
    C: CurveAffine<Base = BF>,
{
    unimplemented!();
}

pub fn load_random_point<'v, BF, C>(ctx: &mut Context<'_, BF>) -> EcPoint<BF>
where
    BF: PrimeField,
    C: CurveAffineExt<Base = BF>,
{
    unimplemented!();
}

// need to supply an extra generic `C` implementing `CurveAffine` trait in order
// to generate random witness points on the curve in question Using Simultaneou, FC::FieldPoints 2^w-Ary Method, see https://www.bmoeller.de/pdf/multiexp-sac2001.pdf
// Random Accumlation point trick learned from halo2wrong: https://hackmd.io/ncuKqRXzR-Cw-Au2fGzsMg?view
// Input:
// - `scalars` is vector of same length as `P`
// - each `scalar` in `scalars` satisfies same assumptions as in
//   `scalar_multiply` above
pub fn multi_scalar_multiply<'v, BF: PrimeField, C>(
    ctx: &mut Context<'_, BF>,
    P: &[EcPoint<BF>],
    scalars: &[<ScalarFieldChip<C> as FieldChip<BF>>::FieldPoint],
    max_bits: usize,
    window_bits: usize,
) -> EcPoint<BF>
where
    C: CurveAffineExt<Base = BF>,
{
    unimplemented!();
}

pub fn get_naf(mut exp: Vec<u64>) -> Vec<i8> {
    // https://en.wikipedia.org/wiki/Non-adjacent_form
    // NAF for exp:
    let mut naf: Vec<i8> = Vec::with_capacity(64 * exp.len());
    let len = exp.len();

    // generate the NAF for exp
    for idx in 0..len {
        let mut e: u64 = exp[idx];
        for _ in 0..64 {
            if e & 1 == 1 {
                let z = 2i8 - (e % 4) as i8;
                e /= 2;
                if z == -1 {
                    e += 1;
                }
                naf.push(z);
            } else {
                naf.push(0);
                e /= 2;
            }
        }
        if e != 0 {
            assert_eq!(e, 1);
            let mut j = idx + 1;
            while j < exp.len() && exp[j] == u64::MAX {
                exp[j] = 0;
                j += 1;
            }
            if j < exp.len() {
                exp[j] += 1;
            } else {
                exp.push(1);
            }
        }
    }
    if exp.len() != len {
        assert_eq!(len, exp.len() + 1);
        assert!(exp[len] == 1);
        naf.push(1);
    }
    naf
}

#[derive(Clone, Debug)]
pub struct EccChip<BF: PrimeField> {
    gate: FlexGateConfig<BF>,
}

impl<BF: PrimeField> EccChip<BF> {
    pub fn construct(gate: FlexGateConfig<BF>) -> Self {
        Self { gate }
    }

    pub fn load_private(
        &self,
        ctx: &mut Context<'_, BF>,
        point: (Value<BF>, Value<BF>),
    ) -> EcPoint<BF> {
        unimplemented!()
    }

    /// Does not constrain witness to lie on curve
    pub fn assign_point<'v, C>(&self, ctx: &mut Context<'_, BF>, g: Value<C>) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        unimplemented!()
    }

    pub fn assign_constant_point<'v, C>(&self, ctx: &mut Context<'_, BF>, g: C) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        unimplemented!();
    }

    pub fn load_random_point<'v, C>(&self, ctx: &mut Context<'_, BF>) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        unimplemented!();
    }

    pub fn assert_is_on_curve<'v, C>(&self, ctx: &mut Context<'_, BF>, P: &EcPoint<BF>)
    where
        C: CurveAffine<Base = BF>,
    {
        is_on_curve::<BF, C>(ctx, P)
    }

    pub fn is_on_curve_or_infinity<'v, C>(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &EcPoint<BF>,
    ) -> AssignedValue<BF>
    where
        C: CurveAffine<Base = BF>,
        C::Base: crate::halo2_proofs::ff::PrimeField,
    {
        unimplemented!();
    }

    pub fn negate(&self, ctx: &mut Context<'_, BF>, P: &EcPoint<BF>) -> EcPoint<BF> {
        unimplemented!();
    }

    /// Assumes that P.x != Q.x
    /// If `is_strict == true`, then actually constrains that `P.x != Q.x`
    pub fn add_unequal(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &EcPoint<BF>,
        Q: &EcPoint<BF>,
        is_strict: bool,
    ) -> EcPoint<BF> {
        ec_add_unequal(ctx, P, Q, is_strict)
    }

    /// Assumes that P.x != Q.x
    /// Otherwise will panic
    pub fn sub_unequal(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &EcPoint<BF>,
        Q: &EcPoint<BF>,
        is_strict: bool,
    ) -> EcPoint<BF> {
        ec_sub_unequal(ctx, P, Q, is_strict)
    }

    pub fn double(&self, ctx: &mut Context<'_, BF>, P: &EcPoint<BF>) -> EcPoint<BF> {
        ec_double(ctx, P)
    }

    pub fn is_equal(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &EcPoint<BF>,
        Q: &EcPoint<BF>,
    ) -> AssignedValue<BF> {
        unimplemented!();
    }

    pub fn assert_equal(&self, ctx: &mut Context<'_, BF>, P: &EcPoint<BF>, Q: &EcPoint<BF>) {
        unimplemented!();
    }

    pub fn sum<'b, 'v: 'b, C>(
        &self,
        ctx: &mut Context<'_, BF>,
        points: impl Iterator<Item = &'b EcPoint<BF>>,
    ) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        unimplemented!();
    }
}

impl<BF: PrimeField> EccChip<BF> {
    pub fn select(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &EcPoint<BF>,
        Q: &EcPoint<BF>,
        condition: &AssignedValue<BF>,
    ) -> EcPoint<BF> {
        ec_select(ctx, P, Q, condition)
    }

    pub fn scalar_mult<'v, C>(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &EcPoint<BF>,
        scalar: &<ScalarFieldChip<C> as FieldChip<BF>>::FieldPoint,
        max_bits: usize,
        window_bits: usize,
    ) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
        C::Base: crate::halo2_proofs::ff::PrimeField,
    {
        scalar_multiply::<BF>(ctx, P, scalar, max_bits, window_bits)
    }

    // TODO: put a check in place that scalar is < modulus of C::Scalar
    pub fn variable_base_msm<'v, C>(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &[EcPoint<BF>],
        scalars: &[<ScalarFieldChip<C> as FieldChip<BF>>::FieldPoint],
        max_bits: usize,
        window_bits: usize,
    ) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
        C::Base: crate::halo2_proofs::ff::PrimeField,
    {
        #[cfg(feature = "display")]
        println!("computing length {} MSM", P.len());

        if P.len() <= 25 {
            multi_scalar_multiply::<BF, C>(ctx, P, scalars, max_bits, window_bits)
        } else {
            // let mut radix = (f64::from((max_bits * scalars[0].len()) as u32)
            // / f64::from(P.len() as u32))
            // .sqrt()
            // .floor() as usize;
            // if radix == 0 {
            // radix = 1;
            // }
            let radix = 1;
            pippenger::multi_exp::<BF, C>(ctx, P, scalars, max_bits, radix, window_bits)
        }
    }
}

impl<BF: PrimeField> EccChip<BF> {
    // TODO: put a check in place that scalar is < modulus of C::Scalar
    pub fn fixed_base_scalar_mult<'v, C>(
        &self,
        ctx: &mut Context<'_, BF>,
        point: &C,
        scalar: &<ScalarFieldChip<C> as FieldChip<BF>>::FieldPoint,
        max_bits: usize,
        window_bits: usize,
    ) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        fixed_base::scalar_multiply::<BF, _>(ctx, point, scalar, max_bits, window_bits)
    }

    /// `radix = 0` means auto-calculate
    ///
    /// `clump_factor = 0` means auto-calculate
    ///
    /// The user should filter out base points that are identity beforehand; we
    /// do not separately do this here
    pub fn fixed_base_msm<'v, C>(
        &self,
        ctx: &mut Context<'_, BF>,
        points: &[C],
        scalars: &[Vec<AssignedValue<BF>>],
        max_scalar_bits_per_cell: usize,
        _radix: usize,
        clump_factor: usize,
    ) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        assert_eq!(points.len(), scalars.len());
        #[cfg(feature = "display")]
        println!("computing length {} fixed base msm", points.len());

        fixed_base::msm(
            self,
            ctx,
            points,
            scalars,
            max_scalar_bits_per_cell,
            clump_factor,
        )

        // Empirically does not seem like pippenger is any better for fixed base
        // msm right now, because of the cost of `select_by_indicator`
        // Cell usage becomes around comparable when `points.len() > 100`, and
        // `clump_factor` should always be 4
        // let radix = if radix == 0 {
        // auto calculate
        // (f64::from(FC::FieldType::NUM_BITS) / f64::from(points.len() as
        // u32)).sqrt().ceil() as usize
        // } else {
        // radix
        // };
        // assert!(radix > 0);
        //
        // fixed_base_pippenger::multi_exp::<BF, C>(
        // self.field_chip,
        // ctx,
        // points,
        // scalars,
        // max_scalar_bits_per_cell,
        // radix,
        // clump_factor,
        // )
    }
}

#[cfg(test)]
pub(crate) mod tests;
