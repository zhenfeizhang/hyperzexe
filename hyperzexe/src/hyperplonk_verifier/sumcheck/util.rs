pub(crate) fn compute_factored_evals<L>(r: &[L::LoadedScalar], righ_num_var: usize) {}
pub(crate) fn eq_poly<L>(r: &[L::LoadedScalar]) {}

pub(crate) fn barycentric_weights<L>(points: &[L::LoadedScalar]) -> Vec<F> {}
pub(crate) fn barycentric_interpolate<L>(
    weights: &[L::LoadedScalar],
    points: &[L::LoadedScalar],
    evals: &[L::LoadedScalar],
    x: &L::LoadedScalar,
) -> F {
}

pub fn evaluate<F: PrimeField>(
    expression: &Expression<L::LoadedScalar>,
    num_vars: usize,
    evals: &HashMap<Query, L::LoadedScalar>,
    challenges: &[L::LoadedScalar],
    ys: &[&[L::LoadedScalar]],
    x: &[L::LoadedScalar],
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
