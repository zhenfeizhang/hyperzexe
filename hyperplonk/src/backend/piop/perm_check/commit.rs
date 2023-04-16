use std::sync::Arc;

use halo2_curves::group::ff::{BatchInvert, PrimeField};
use itertools::Itertools;

use crate::{
    backend::{
        poly::multilinear::MultilinearPolynomial,
        util::{arithmetic::steps_by, end_timer, parallel::parallelize, start_timer},
    },
    Error,
};

/// Compute multilinear fractional polynomial s.t. frac(x) = f1(x) * ... * fk(x)
/// / (g1(x) * ... * gk(x)) for all x \in {0,1}^n
/// - input:
///     - max_degree: maximum degree of the constraints
///     - permutation_polys: [(poly_idx, next_for_poly_idx(x)), ...]
///     - polys: all preprocessed, instance and witness polynomials
pub(super) fn compute_frac_poly<F: PrimeField>(
    max_degree: usize,
    permutation_polys: &[(usize, Arc<MultilinearPolynomial<F>>)],
    polys: &[Arc<MultilinearPolynomial<F>>],
    beta: &F,
    gamma: &F,
) -> Result<Vec<Arc<MultilinearPolynomial<F>>>, Error> {
    let start = start_timer(|| "compute frac(x)");

    let chunk_size = max_degree - 1;
    let polys = polys.into_iter().collect_vec();
    let permutation_polys = permutation_polys.into_iter().collect_vec();
    let num_vars = polys[0].num_vars();

    let products = permutation_polys
        .chunks(chunk_size)
        .enumerate()
        .map(|(chunk_idx, permutation_polys)| {
            let mut product = vec![F::one(); 1 << num_vars];

            for (poly, permutation_poly) in permutation_polys.iter() {
                parallelize(&mut product, |(product, start)| {
                    for ((product, value), permutation) in product
                        .iter_mut()
                        .zip(polys[*poly][start..].iter())
                        .zip(permutation_poly[start..].iter())
                    {
                        *product *= (*beta * permutation) + gamma + value;
                    }
                });
            }

            product.iter_mut().batch_invert();

            for ((poly, _), idx) in permutation_polys.iter().zip(chunk_idx * chunk_size..) {
                let id_offset = idx << num_vars;
                parallelize(&mut product, |(product, start)| {
                    for ((product, value), beta_id) in product
                        .iter_mut()
                        .zip(polys[*poly][start..].iter())
                        .zip(steps_by(F::from((id_offset + start) as u64) * beta, *beta))
                    {
                        *product *= beta_id + gamma + value;
                    }
                });
            }

            Arc::new(MultilinearPolynomial::new(product))
        })
        .collect_vec();
    end_timer(start);
    Ok(products)
}

/// Compute the product polynomial `prod(x)` such that
/// `prod(x) = [(1-xn)*frac(0, x1, ..., x(n-1)) + xn*prod(0, x1, ..., x(n-1))] *
/// [(1-xn)*frac(1, x1, ..., x(n-1)) + xn*prod(1, x1, ..., x(n-1))]` on the
/// boolean hypercube {0,1}^n
///
/// The caller needs to check num_vars() matches in f and g
/// Cost: linear in N.
pub(super) fn compute_product_poly<F: PrimeField>(
    frac_polys: &[Arc<MultilinearPolynomial<F>>],
) -> Result<Vec<Arc<MultilinearPolynomial<F>>>, Error> {
    let start = start_timer(|| "compute prod(x)");
    let num_vars = frac_polys[0].num_vars();
    let prod_polys = frac_polys
        .iter()
        .map(|frac_poly| {
            let frac_evals = &frac_poly.evals();

            // ===================================
            // prod(x)
            // ===================================
            //
            // `prod(x)` can be computed via recursing the following formula for 2^n-1
            // times
            //
            // `prod(x_1, ..., x_n) :=
            //      p1 = [(1 - xn) * frac(0, x1, ..., x(n-1)) + xn * prod(0, x1, ..., x(n-1))] *
            //      p2 = [(1 - xn) * frac(1, x1, ..., x(n-1)) + xn * prod(1, x1, ..., x(n-1))]`
            //
            // At any given step, the right hand side of the equation
            // is available via either frac_x or the current view of prod_x
            let mut prod_x_evals = vec![];
            for x in 0..(1 << num_vars) - 1 {
                // sign will decide if the evaluation should be looked up from frac_x or
                // prod_x; x_zero_index is the index for the evaluation (x_2, ..., x_n,
                // 0); x_one_index is the index for the evaluation (x_2, ..., x_n, 1);
                let (x_zero_index, x_one_index, sign) = get_index(x, num_vars);
                if !sign {
                    prod_x_evals.push(frac_evals[x_zero_index] * frac_evals[x_one_index]);
                } else {
                    // sanity check: if we are trying to look up from the prod_x_evals table,
                    // then the target index must already exist
                    if x_zero_index >= prod_x_evals.len() || x_one_index >= prod_x_evals.len() {
                        return Err(Error::InvalidSnark(
                            "target index must already exist".to_string(),
                        ));
                    }
                    prod_x_evals.push(prod_x_evals[x_zero_index] * prod_x_evals[x_one_index]);
                }
            }

            assert_eq!(prod_x_evals[prod_x_evals.len() - 1], F::one());
            // prod(1, 1, ..., 1) := 0
            prod_x_evals.push(F::zero());

            Ok(Arc::new(MultilinearPolynomial::new(prod_x_evals)))
        })
        .collect::<Result<Vec<_>, _>>()?;
    end_timer(start);
    Ok(prod_polys)
}

/// Compute p1(x) and p2(x) such that
/// p1(x_1, .., x_n) = (1 - x_n) * frac(0, x_1, ..., x_(n-1))
///                  + x_n * prod(0, x_1, ..., x_(n - 1))
///                  = sum_y ((1 - x_n) * eq((0, x_1, ..., x_(n-1)), y) * frac(0, x_1, ..., x_(n - 1))
///                              + x_n * eq((0, x_1, ..., x_(n-1)), y) * prod(0, x_1, ..., x_(n - 1))
/// p2(x_1, .., x_n) = (1 - x_n) * frac(1, x_1, ..., x_(n-1))
///                  + x_n * prod(1, x_1, ..., x_(n - 1))
///                  = sum_y ((1 - x_n) * eq((1, x_1, ..., x_(n-1)), y) * frac(1, x_1, ..., x_(n - 1))
///                             + x_n * eq((1, x_1, ..., x_(n-1)), y) * prod(1, x_1, ..., x_(n - 1))
/// 
pub(super) fn compute_p1_p2_poly<F: PrimeField>(
    frac_poly: &MultilinearPolynomial<F>,
    prod_poly: &MultilinearPolynomial<F>,
) -> (Arc<MultilinearPolynomial<F>>, Arc<MultilinearPolynomial<F>>) {
    let start = start_timer(|| "compute p1(x) and p2(x)");
    let num_vars = frac_poly.num_vars();
    let mut p1_evals = vec![F::zero(); 1 << num_vars];
    let mut p2_evals = vec![F::zero(); 1 << num_vars];
    for x in 0..1 << num_vars {
        let (x0, x1, sign) = get_index(x, num_vars);
        if !sign {
            p1_evals[x] = frac_poly.evals()[x0];
            p2_evals[x] = frac_poly.evals()[x1];
        } else {
            p1_evals[x] = prod_poly.evals()[x0];
            p2_evals[x] = prod_poly.evals()[x1];
        }
    }
    let p1 = Arc::new(MultilinearPolynomial::new(p1_evals));
    let p2 = Arc::new(MultilinearPolynomial::new(p2_evals));
    end_timer(start);
    (p1, p2)
}

// Input index
// - `i := (i_{n-1}, ..., i_0)`,
// - `num_vars := n`
// return three elements:
// - `x0 := (i_{n-2}, ..., i_0, 0)`
// - `x1 := (i_{n-2}, ..., i_0, 1)`
// - `sign := i_{n-1}`
#[inline]
pub fn get_index(i: usize, num_vars: usize) -> (usize, usize, bool) {
    let x0 = (i << 1) & ((1 << num_vars) - 1);
    let x1 = x0 | 1;
    let sign = (i >> (num_vars - 1)) & 1 == 1;

    (x0, x1, sign)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::{
        poly::multilinear::{arc_mpoly, arc_mpoly_fr, MultilinearPolynomial},
        util::arithmetic::{product, sum},
    };
    use halo2_curves::{bn256::Fr as F, group::ff::Field};
    #[test]
    fn test_compute_permutation_help_polys() {
        // Prepare test data
        let max_degree = 3;
        let chunk_size = max_degree - 1;
        let num_vars = 2;
        let beta = F::from(2 as u64); // Provide a valid beta value of type F
        let gamma = F::from(3 as u64); // Provide a valid gamma value of type F

        // Provide appropriate MultilinearPolynomial values
        let polys = vec![arc_mpoly!(F, 0, 1, 2, 3), arc_mpoly!(F, 1, 1, 2, 3)];
        let ids = {
            let mut ids = vec![];
            for i in 0..chunk_size {
                let shift = (i * (1 << num_vars)) as u64;
                let s_id_vec = (shift..shift + (1u64 << num_vars))
                    .map(F::from)
                    .collect_vec();
                ids.push(s_id_vec);
            }
            ids
        };

        let at = |i: usize, j: usize| -> F { ids[i][j] };
        let permutation_polys = vec![
            (0, arc_mpoly_fr![at(0, 0), at(1, 0), at(1, 2), at(1, 3)]),
            (1, arc_mpoly_fr![at(1, 1), at(0, 1), at(0, 2), at(0, 3)]),
        ];

        // Compute `frac`
        let fracs = compute_frac_poly(max_degree, &permutation_polys, &polys, &beta, &gamma)
            .expect("Error occurred during compute_frac_poly");
        let frac_evals = vec![fracs[0].evals().to_vec()];

        macro_rules! fx {
            ($i:expr, $j:expr) => {
                polys[$i].evals()[$j] + beta * at($i % chunk_size, $j) + gamma
            };
        }
        macro_rules! gx {
            ($i:expr, $j:expr) => {
                polys[$i].evals()[$j] + beta * permutation_polys[$i].1.evals()[$j] + gamma
            };
        }

        // Define expected `frac`
        let want_fracs = vec![vec![
            fx!(0, 0) * fx!(1, 0) * (gx!(0, 0) * gx!(1, 0)).invert().unwrap(),
            fx!(0, 1) * fx!(1, 1) * (gx!(0, 1) * gx!(1, 1)).invert().unwrap(),
            fx!(0, 2) * fx!(1, 2) * (gx!(0, 2) * gx!(1, 2)).invert().unwrap(),
            fx!(0, 3) * fx!(1, 3) * (gx!(0, 3) * gx!(1, 3)).invert().unwrap(),
        ]];

        // Check if `frac` is as expected

        assert_eq!(frac_evals[0], want_fracs[0]);
        assert_eq!(product(frac_evals[0].clone()), F::one());

        // Compute `prod`
        let prods =
            compute_product_poly(&fracs).expect("Error occurred during compute_product_poly");
        let prod_evals = vec![prods[0].evals().to_vec()];

        // Check if `prod` is as expected
        let want_prods = vec![vec![
            frac_evals[0][0] * frac_evals[0][1],
            frac_evals[0][2] * frac_evals[0][3],
            frac_evals[0][0] * frac_evals[0][1] * frac_evals[0][2] * frac_evals[0][3],
        ]];
        assert_eq!(prod_evals[0][..3], want_prods[0]);

        // Compute `p1` and `p2`
        let mut p1_evals = vec![];
        let mut p2_evals = vec![];
        let mut p1s = vec![];
        let mut p2s = vec![];
        for (frac, prod) in fracs.iter().zip(prods.iter()) {
            let (p1, p2) = compute_p1_p2_poly(frac, prod);
            let p1 = Arc::try_unwrap(p1).unwrap();
            let p2 = Arc::try_unwrap(p2).unwrap();
            p1_evals.push(p1.evals().to_vec());
            p2_evals.push(p2.evals().to_vec());
            p1s.push(p1);
            p2s.push(p2);
        }
        let want_p1s = vec![vec![
            frac_evals[0][0],
            frac_evals[0][2],
            prod_evals[0][0],
            prod_evals[0][2],
        ]];
        assert_eq!(p1_evals[0], want_p1s[0]);

        let want_p2s = vec![vec![
            frac_evals[0][1],
            frac_evals[0][3],
            prod_evals[0][1],
            prod_evals[0][3],
        ]];
        assert_eq!(p2_evals[0], want_p2s[0]);

        // Test expressions
        let point = vec![F::from(97), F::from(41)];
        let frac_left = MultilinearPolynomial::eq_xy_shifted_ab(&point, false, false).into_evals();
        let frac_right = MultilinearPolynomial::eq_xy_shifted_ab(&point, false, true).into_evals();
        let prod_left = MultilinearPolynomial::eq_xy_shifted_ab(&point, true, false).into_evals();
        let prod_right = MultilinearPolynomial::eq_xy_shifted_ab(&point, true, true).into_evals();
        let want_p1_eval = p1s[0].evaluate(&point);
        let want_p2_eval = p2s[0].evaluate(&point);
        let p1_eval =
            sum((0..4).map(|i| frac_left[i] * frac_evals[0][i] + prod_left[i] * prod_evals[0][i]));
        let p2_eval = sum(
            (0..4).map(|i| frac_right[i] * frac_evals[0][i] + prod_right[i] * prod_evals[0][i])
        );
        assert_eq!(want_p1_eval, p1_eval);
        assert_eq!(want_p2_eval, p2_eval);
    }
}
