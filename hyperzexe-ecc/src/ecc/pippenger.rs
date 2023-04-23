use super::{
    ec_add_unequal, ec_double, ec_select, ec_select_from_bits, ec_sub_unequal, load_random_point,
    EcPoint,
};
use crate::fields::{FieldChip, Selectable};
use halo2_base::{
    gates::GateInstructions,
    utils::{CurveAffineExt, PrimeField},
    AssignedValue, Context,
};

// Reference: https://jbootle.github.io/Misc/pippenger.pdf

// Reduction to multi-products
// Output:
// * new_points: length `points.len() * radix`
// * new_bool_scalars: 2d array `ceil(scalar_bits / radix)` by `points.len() *
//   radix`
pub fn decompose<'v, F>(
    ctx: &mut Context<'_, F>,
    points: &[EcPoint<F>],
    scalars: &[Vec<AssignedValue<F>>],
    max_scalar_bits_per_cell: usize,
    radix: usize,
) -> (Vec<EcPoint<F>>, Vec<Vec<AssignedValue<F>>>)
where
    F: PrimeField,
{
    unimplemented!();
}

// Given points[i] and bool_scalars[j][i],
// compute G'[j] = sum_{i=0..points.len()} points[i] * bool_scalars[j][i]
// output is [ G'[j] + rand_point ]_{j=0..bool_scalars.len()}, rand_point
pub fn multi_product<'v, F: PrimeField, C>(
    ctx: &mut Context<'_, F>,
    points: &[EcPoint<F>],
    bool_scalars: &[Vec<AssignedValue<F>>],
    clumping_factor: usize,
) -> (Vec<EcPoint<F>>, EcPoint<F>)
where
    C: CurveAffineExt<Base = F>,
{
    unimplemented!();
}

pub fn multi_exp<'v, F: PrimeField, C>(
    ctx: &mut Context<'_, F>,
    points: &[EcPoint<F>],
    scalars: &[Vec<AssignedValue<F>>],
    max_scalar_bits_per_cell: usize,
    radix: usize,
    clump_factor: usize,
) -> EcPoint<F>
where
    C: CurveAffineExt<Base = F>,
{
    unimplemented!();
}
