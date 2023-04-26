use std::{collections::HashMap, iter};

use halo2_proofs::curves::CurveAffine;
use hyperplonk::backend::util::expression::{CommonPolynomial, Expression, Query};

use crate::halo2_verifier::loader::Loader;

pub fn eq_eval<C: CurveAffine, L: Loader<C>>(r: &[L::LoadedScalar]) -> Vec<L::LoadedScalar> {}

pub fn compute_factored_evals<C: CurveAffine, L: Loader<C>>(
    r: &[L::LoadedScalar],
    right_num_vars: usize,
) -> (Vec<L::LoadedScalar>, Vec<L::LoadedScalar>) {
    let ell = r.len();

    let L = eq_eval(&r[right_num_vars..ell]);
    let R = eq_eval(&r[..right_num_vars]);

    (L, R)
}

pub fn barycentric_weights<C: CurveAffine, L: Loader<C>>(
    points: &[L::LoadedScalar],
) -> Vec<L::LoadedScalar> {
    unimplemented!();
}

pub fn barycentric_interpolate<C: CurveAffine, L: Loader<C>>(
    weights: &[L::LoadedScalar],
    points: &[L::LoadedScalar],
    evals: &[L::LoadedScalar],
    x: &L::LoadedScalar,
) -> L::LoadedScalar {
    unimplemented!();
}

pub fn evaluate<C: CurveAffine, L: Loader<C>>(
    expression: &Expression<L::LoadedScalar>,
    num_vars: usize,
    evals: &HashMap<Query, L::LoadedScalar>,
    challenges: &[L::LoadedScalar],
    ys: &[&[L::LoadedScalar]],
    x: &[L::LoadedScalar],
) -> L::LoadedScalar {
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
            CommonPolynomial::Identity(idx) => {
                let loader = x[0].loader();
                loader.load_constant(C::Scalar::from((idx << num_vars) as u64)) + identity
            },
        },
        &|query| evals[&query],
        &|idx| challenges[idx],
        &|scalar| -scalar,
        &|lhs, rhs| lhs + &rhs,
        &|lhs, rhs| lhs * &rhs,
        &|value, scalar| scalar * value,
    )
}

pub fn lagrange_eval<C: CurveAffine, L: Loader<C>>(
    x: &[L::LoadedScalar],
    b: usize,
) -> L::LoadedScalar {
    assert!(!x.is_empty());

    let loader = x[0].loader();
    let one = loader.load_one();
    loader.product(x.iter().enumerate().map(
        |(idx, x_i)| {
            if b.nth_bit(idx) {
                *x_i
            } else {
                one - x_i
            }
        },
    ))
}

pub fn eq_xy_eval<C: CurveAffine, L: Loader<C>>(
    x: &[L::LoadedScalar],
    y: &[L::LoadedScalar],
) -> L::LoadedScalar {
    assert!(!x.is_empty());
    assert_eq!(x.len(), y.len());

    let loader = x[0].loader();
    let one = loader.load_one();
    loader.product(
        x.iter()
            .zip(y)
            .map(|(x_i, y_i)| (*x_i * y_i).double() + one - x_i - y_i),
    )
}

/// Compute ((1 - yn) * (1 - a) + yn * a ) eq ((b, y1, ..., yn), x)
pub fn eq_xy_eval_shift_ab<C: CurveAffine, L: Loader<C>>(
    x: &[L::LoadedScalar],
    y: &[L::LoadedScalar],
    a: bool,
    b: bool,
) -> L::LoadedScalar {
    assert!(!x.is_empty());
    let loader = x[0].loader();
    let zero = loader.load_zero();
    let one = loader.load_one();
    let y_last = if a {
        y[y.len() - 1]
    } else {
        one - y[y.len() - 1]
    };
    let y_0 = if b { one } else { zero };
    let y = iter::once(y_0)
        .chain(y[..y.len() - 1].to_vec())
        .collect_vec();
    y_last * eq_xy_eval(x, &y)
}

fn identity_eval<C: CurveAffine, L: Loader<C>>(x: &[L::LoadedScalar]) -> C::Scalar {
    let two_power = C::Scalar::from(2u64);
    let two_powers = Vec::with_capacity(x.len());
    for i in 0..x.len() {
        two_powers.push(two_power);
        two_power = two_power.double();
    }
    let loader = x[0].loader();
    loader.inner_product(x, &two_powers)
}
