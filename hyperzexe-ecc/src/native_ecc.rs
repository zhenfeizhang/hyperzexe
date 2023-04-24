#![allow(non_snake_case)]
use crate::halo2_proofs::{
    arithmetic::CurveAffine,
    circuit::Value,
    group::{Curve, Group},
};
use halo2_base::{
    gates::{flex_gate::FlexGateConfig, GateInstructions, RangeInstructions},
    halo2_proofs::group::prime::PrimeCurveAffine,
    utils::{modulus, CurveAffineExt, PrimeField},
    AssignedValue, Context,
    QuantumCell::{Constant, Existing},
};
use halo2_ecc::fields::{fp::FpConfig, FieldChip, PrimeFieldChip, Selectable};
use itertools::Itertools;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::marker::PhantomData;

pub mod fixed_base;
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
pub fn ec_add_unequal<'v, BF: PrimeField, Gate: GateInstructions<BF>>(
    gate: &Gate,
    ctx: &mut Context<'_, BF>,
    P: &EcPoint<BF>,
    Q: &EcPoint<BF>,
    is_strict: bool,
) -> EcPoint<BF> {
    if is_strict {
        // constrains that P.x != Q.x
        let x_is_equal = gate.is_equal(ctx, Existing(&P.x), Existing(&Q.x));
        gate.assert_is_const(ctx, &x_is_equal, BF::ZERO);
    }

    let dx = gate.sub(ctx, Existing(&Q.x), Existing(&P.x));
    let dy = gate.sub(ctx, Existing(&Q.y), Existing(&P.y));
    let lambda = gate.div_unsafe(ctx, Existing(&dy), Existing(&dx));

    //  x_3 = lambda^2 - x_1 - x_2 (mod p)
    let lambda_sq = gate.mul(ctx, Existing(&lambda), Existing(&lambda));
    let lambda_sq_minus_px = gate.sub(ctx, Existing(&lambda_sq), Existing(&P.x));
    let x_3 = gate.sub(ctx, Existing(&lambda_sq_minus_px), Existing(&Q.x));

    //  y_3 = lambda (x_1 - x_3) - y_1 mod p
    let dx_13 = gate.sub(ctx, Existing(&P.x), Existing(&x_3));
    let lambda_dx_13 = gate.mul(ctx, Existing(&lambda), Existing(&dx_13));
    let y_3 = gate.sub(ctx, Existing(&lambda_dx_13), Existing(&P.y));

    EcPoint::construct(x_3, y_3)
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
pub fn ec_sub_unequal<'v, BF: PrimeField, Gate: GateInstructions<BF>>(
    gate: &Gate,
    ctx: &mut Context<'_, BF>,
    P: &EcPoint<BF>,
    Q: &EcPoint<BF>,
    is_strict: bool,
) -> EcPoint<BF> {
    if is_strict {
        let x_is_equal = gate.is_equal(ctx, Existing(&P.x), Existing(&Q.x));
        gate.assert_is_const(ctx, &x_is_equal, BF::ZERO);
    }

    let dx = gate.sub(ctx, Existing(&Q.x), Existing(&P.x));
    let dy = gate.add(ctx, Existing(&Q.y), Existing(&P.y));

    let neg_dy = gate.neg(ctx, Existing(&dy));
    let lambda = gate.div_unsafe(ctx, Existing(&neg_dy), Existing(&dx));

    // (x_2 - x_1) * lambda + y_2 + y_1 = 0 (mod p)
    let lambda_dx = gate.mul(ctx, Existing(&lambda), Existing(&dx));
    let lambda_dx_plus_dy = gate.add(ctx, Existing(&lambda_dx), Existing(&dy));
    gate.assert_is_const(ctx, &lambda_dx_plus_dy, BF::ZERO);

    //  x_3 = lambda^2 - x_1 - x_2 (mod p)
    let lambda_sq = gate.mul(ctx, Existing(&lambda), Existing(&lambda));
    let lambda_sq_minus_px = gate.sub(ctx, Existing(&lambda_sq), Existing(&P.x));
    let x_3 = gate.sub(ctx, Existing(&lambda_sq_minus_px), Existing(&Q.x));

    //  y_3 = lambda (x_1 - x_3) - y_1 mod p
    let dx_13 = gate.sub(ctx, Existing(&P.x), Existing(&x_3));
    let lambda_dx_13 = gate.mul(ctx, Existing(&lambda), Existing(&dx_13));
    let y_3 = gate.sub(ctx, Existing(&lambda_dx_13), Existing(&P.y));

    EcPoint::construct(x_3, y_3)
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
pub fn ec_double<'v, BF: PrimeField, Gate: GateInstructions<BF>>(
    gate: &Gate,
    ctx: &mut Context<'_, BF>,
    P: &EcPoint<BF>,
) -> EcPoint<BF> {
    // removed optimization that computes `2 * lambda` while assigning witness to
    // `lambda` simultaneously, in favor of readability. The difference is just
    // copying `lambda` once
    let two = BF::from(2u64);
    let two_y = gate.mul(ctx, Existing(&P.y), Constant(two.clone()));
    let three = BF::from(3u64);
    let three_x = gate.mul(ctx, Existing(&P.x), Constant(three));
    let three_x_sq = gate.mul(ctx, Existing(&three_x), Existing(&P.x));
    let lambda = gate.div_unsafe(ctx, Existing(&three_x_sq), Existing(&two_y));

    // x_3 = lambda^2 - 2 x % p
    let lambda_sq = gate.mul(ctx, Existing(&lambda), Existing(&lambda));
    let two_x = gate.mul(ctx, Existing(&P.x), Constant(two.clone()));
    let x_3 = gate.sub(ctx, Existing(&lambda_sq), Existing(&two_x));

    // y_3 = lambda (x - x_3) - y % p
    let dx = gate.sub(ctx, Existing(&P.x), Existing(&x_3));
    let lambda_dx = gate.mul(ctx, Existing(&lambda), Existing(&dx));
    let y_3 = gate.sub(ctx, Existing(&lambda_dx), Existing(&P.y));

    EcPoint::construct(x_3, y_3)
}

pub fn ec_select<'v, BF: PrimeField, Gate: GateInstructions<BF>>(
    gate: &Gate,
    ctx: &mut Context<'_, BF>,
    P: &EcPoint<BF>,
    Q: &EcPoint<BF>,
    sel: &AssignedValue<BF>,
) -> EcPoint<BF> {
    let Rx = gate.select(ctx, Existing(&P.x), Existing(&Q.x), Existing(sel));
    let Ry = gate.select(ctx, Existing(&P.y), Existing(&Q.y), Existing(sel));
    EcPoint::construct(Rx, Ry)
}

// takes the dot product of points with sel, where each is intepreted as
// a _vector_
pub fn ec_select_by_indicator<'v, BF: PrimeField, Gate: GateInstructions<BF>>(
    gate: &Gate,
    ctx: &mut Context<'_, BF>,
    points: &[EcPoint<BF>],
    coeffs: &[AssignedValue<BF>],
) -> EcPoint<BF> {
    let x_coords = points.iter().map(|P| Existing(&P.x)).collect::<Vec<_>>();
    let y_coords = points.iter().map(|P| Existing(&P.y)).collect::<Vec<_>>();
    let Rx = gate.select_by_indicator(ctx, x_coords, coeffs);
    let Ry = gate.select_by_indicator(ctx, y_coords, coeffs);
    EcPoint::construct(Rx, Ry)
}

// `sel` is little-endian binary
pub fn ec_select_from_bits<'v, BF: PrimeField, Gate: GateInstructions<BF>>(
    gate: &Gate,
    ctx: &mut Context<'_, BF>,
    points: &[EcPoint<BF>],
    sel: &[AssignedValue<BF>],
) -> EcPoint<BF> {
    let w = sel.len();
    let num_points = points.len();
    assert_eq!(1 << w, num_points);
    let coeffs = gate.bits_to_indicator(ctx, sel);
    ec_select_by_indicator(gate, ctx, points, &coeffs)
}

// computes [scalar] * P on y^2 = x^3 + b
//   * P has order given by the scalar field modulus
pub fn scalar_multiply<'v, BF: PrimeField, SFPoint>(
    ctx: &mut Context<'_, BF>,
    P: &EcPoint<BF>,
    scalar: &SFPoint,
    max_bits: usize,
    window_bits: usize,
) -> EcPoint<BF> {
    unimplemented!();
}

pub fn is_on_curve<'v, C, Gate>(gate: &Gate, ctx: &mut Context<'_, C::Base>, P: &EcPoint<C::Base>)
where
    C: CurveAffine,
    C::Base: PrimeField,
    Gate: GateInstructions<C::Base>,
{
    let lhs = gate.mul(ctx, Existing(&P.y), Existing(&P.y));
    let mut rhs = gate.mul(ctx, Existing(&P.x), Existing(&P.x));
    rhs = gate.mul(ctx, Existing(&rhs), Existing(&P.x));

    let b = gate.load_constant(ctx, C::b());
    rhs = gate.add(ctx, Existing(&rhs), Existing(&b));
    let diff = gate.sub(ctx, Existing(&lhs), Existing(&rhs));
    gate.is_zero(ctx, &diff);
}

pub fn load_random_point<'v, C, Gate>(
    gate: &Gate,
    ctx: &mut Context<'_, C::Base>,
) -> EcPoint<C::Base>
where
    C: CurveAffineExt,
    C::Base: PrimeField,
    Gate: GateInstructions<C::Base>,
{
    let base_point: C = C::CurveExt::random(ChaCha20Rng::from_entropy()).to_affine();
    let (x, y) = base_point.into_coordinates();
    let x_assigned = gate.load_constant(ctx, x);
    let y_assigned = gate.load_constant(ctx, y);
    let base = EcPoint::construct(x_assigned, y_assigned);
    // for above reason we still need to constrain that the witness is on the curve
    is_on_curve::<C, Gate>(gate, ctx, &base);
    base
}

// need to supply an extra generic `C` implementing `CurveAffine` trait in order
// to generate random witness points on the curve in question Using Simultaneou, FC::FieldPoints 2^w-Ary Method, see https://www.bmoeller.de/pdf/multiexp-sac2001.pdf
// Random Accumlation point trick learned from halo2wrong: https://hackmd.io/ncuKqRXzR-Cw-Au2fGzsMg?view
// Input:
// - `scalars` is vector of same length as `P`
// - each `scalar` in `scalars` satisfies same assumptions as in
//   `scalar_multiply` above
pub fn multi_scalar_multiply<'v, BF: PrimeField, SFC>(
    ctx: &mut Context<'_, BF>,
    P: &[EcPoint<BF>],
    scalars: &[SFC::FieldPoint],
    max_bits: usize,
    window_bits: usize,
) -> EcPoint<BF>
where
    SFC: FieldChip<BF>,
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
pub struct EccChip<BF: PrimeField, Gate: GateInstructions<BF>> {
    pub gate: Gate,
    _marker: PhantomData<BF>,
}

impl<BF: PrimeField, Gate: GateInstructions<BF>> EccChip<BF, Gate> {
    pub fn construct(gate: Gate) -> Self {
        Self {
            gate,
            _marker: PhantomData,
        }
    }

    pub fn load_private(
        &self,
        ctx: &mut Context<'_, BF>,
        point: (Value<BF>, Value<BF>),
    ) -> EcPoint<BF> {
        let x_assigned = self.gate.load_witness(ctx, point.0);
        let y_assigned = self.gate.load_witness(ctx, point.1);
        EcPoint::construct(x_assigned, y_assigned)
    }

    /// Does not constrain witness to lie on curve
    pub fn assign_point<'v, C>(&self, ctx: &mut Context<'_, BF>, g: Value<C>) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        let (x, y) = g.map(|g| g.into_coordinates()).unzip();
        let x_assigned = self.gate.load_witness(ctx, x);
        let y_assigned = self.gate.load_witness(ctx, y);
        EcPoint::construct(x_assigned, y_assigned)
    }

    pub fn assign_constant_point<'v, C>(&self, ctx: &mut Context<'_, BF>, g: C) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        let (x, y) = g.into_coordinates();
        let x_assigned = self.gate.load_constant(ctx, x);
        let y_assigned = self.gate.load_constant(ctx, y);
        EcPoint::construct(x_assigned, y_assigned)
    }

    pub fn load_random_point<'v, C>(&self, ctx: &mut Context<'_, BF>) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        load_random_point::<C, Gate>(&self.gate, ctx)
    }

    pub fn assert_is_on_curve<'v, C>(&self, ctx: &mut Context<'_, BF>, P: &EcPoint<BF>)
    where
        C: CurveAffine<Base = BF>,
    {
        is_on_curve::<C, Gate>(&self.gate, ctx, P)
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
        let lhs = self.gate.mul(ctx, Existing(&P.y), Existing(&P.y));
        let mut rhs = self.gate.mul(ctx, Existing(&P.x), Existing(&P.x));
        rhs = self.gate.mul(ctx, Existing(&rhs), Existing(&P.x));

        let b = self.gate.load_constant(ctx, C::b());
        rhs = self.gate.add(ctx, Existing(&rhs), Existing(&b));
        let diff = self.gate.sub(ctx, Existing(&lhs), Existing(&rhs));

        let is_on_curve = self.gate.is_zero(ctx, &diff);
        let x_is_zero = self.gate.is_zero(ctx, &P.x);
        let y_is_zero = self.gate.is_zero(ctx, &P.y);

        self.gate.or_and(
            ctx,
            Existing(&is_on_curve),
            Existing(&x_is_zero),
            Existing(&y_is_zero),
        )
    }

    pub fn negate(&self, ctx: &mut Context<'_, BF>, P: &EcPoint<BF>) -> EcPoint<BF> {
        EcPoint::construct(P.x.clone(), self.gate.neg(ctx, Existing(&P.y)))
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
        ec_add_unequal(&self.gate, ctx, P, Q, is_strict)
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
        ec_sub_unequal(&self.gate, ctx, P, Q, is_strict)
    }

    pub fn double(&self, ctx: &mut Context<'_, BF>, P: &EcPoint<BF>) -> EcPoint<BF> {
        ec_double(&self.gate, ctx, P)
    }

    pub fn is_equal(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &EcPoint<BF>,
        Q: &EcPoint<BF>,
    ) -> AssignedValue<BF> {
        let x_is_equal = self.gate.is_equal(ctx, Existing(&P.x), Existing(&Q.x));
        let y_is_equal = self.gate.is_equal(ctx, Existing(&P.y), Existing(&Q.y));
        self.gate
            .and(ctx, Existing(&x_is_equal), Existing(&y_is_equal))
    }

    pub fn assert_equal(&self, ctx: &mut Context<'_, BF>, P: &EcPoint<BF>, Q: &EcPoint<BF>) {
        self.gate.assert_equal(ctx, Existing(&P.x), Existing(&Q.x));
        self.gate.assert_equal(ctx, Existing(&P.y), Existing(&Q.y));
    }

    pub fn sum<'b, 'v: 'b, C>(
        &self,
        ctx: &mut Context<'_, BF>,
        points: impl Iterator<Item = &'b EcPoint<BF>>,
    ) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
    {
        //? Why do we need to add a random point?
        let rand_point = self.load_random_point::<C>(ctx);
        let mut acc = rand_point.clone();
        for point in points {
            acc = self.add_unequal(ctx, &acc, point, true);
        }
        self.sub_unequal(ctx, &acc, &rand_point, true)
    }
}

impl<BF: PrimeField, Gate: GateInstructions<BF>> EccChip<BF, Gate> {
    pub fn select(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &EcPoint<BF>,
        Q: &EcPoint<BF>,
        condition: &AssignedValue<BF>,
    ) -> EcPoint<BF> {
        ec_select(&self.gate, ctx, P, Q, condition)
    }

    pub fn scalar_mult<'v, SFPoint>(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &EcPoint<BF>,
        scalar: &SFPoint,
        max_bits: usize,
        window_bits: usize,
    ) -> EcPoint<BF> {
        scalar_multiply::<BF, SFPoint>(ctx, P, scalar, max_bits, window_bits)
    }

    // TODO: put a check in place that scalar is < modulus of C::Scalar
    pub fn variable_base_msm<'v, C, SFC>(
        &self,
        ctx: &mut Context<'_, BF>,
        P: &[EcPoint<BF>],
        scalars: &[SFC::FieldPoint],
        max_bits: usize,
        window_bits: usize,
    ) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
        SFC: FieldChip<BF, FieldType = C::Scalar>,
    {
        #[cfg(feature = "display")]
        println!("computing length {} MSM", P.len());

        if P.len() <= 25 {
            multi_scalar_multiply::<BF, SFC>(ctx, P, scalars, max_bits, window_bits)
        } else {
            // let mut radix = (f64::from((max_bits * scalars[0].len()) as u32)
            // / f64::from(P.len() as u32))
            // .sqrt()
            // .floor() as usize;
            // if radix == 0 {
            // radix = 1;
            // }
            let radix = 1;
            pippenger::multi_exp::<BF, C, SFC>(ctx, P, scalars, max_bits, radix, window_bits)
        }
    }
}

impl<BF: PrimeField, Gate: GateInstructions<BF>> EccChip<BF, Gate> {
    // TODO: put a check in place that scalar is < modulus of C::Scalar
    pub fn fixed_base_scalar_mult<'v, C, SFC>(
        &self,
        ctx: &mut Context<'_, BF>,
        point: &C,
        scalar: &SFC::FieldPoint,
        max_bits: usize,
        window_bits: usize,
    ) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
        SFC: FieldChip<BF, FieldType = C::Scalar>,
    {
        fixed_base::scalar_multiply::<C, SFC>(ctx, point, scalar, max_bits, window_bits)
    }

    /// `radix = 0` means auto-calculate
    ///
    /// `clump_factor = 0` means auto-calculate
    ///
    /// The user should filter out base points that are identity beforehand; we
    /// do not separately do this here
    pub fn fixed_base_msm<'v, C, SFC>(
        &self,
        ctx: &mut Context<'_, BF>,
        points: &[C],
        scalars: &[SFC::FieldPoint],
        max_scalar_bits_per_cell: usize,
        _radix: usize,
        clump_factor: usize,
    ) -> EcPoint<BF>
    where
        C: CurveAffineExt<Base = BF>,
        SFC: FieldChip<BF, FieldType = C::Scalar>,
    {
        assert_eq!(points.len(), scalars.len());
        #[cfg(feature = "display")]
        println!("computing length {} fixed base msm", points.len());

        fixed_base::msm::<C, SFC, Gate>(
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
