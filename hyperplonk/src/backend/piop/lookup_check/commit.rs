use ark_std::{end_timer, start_timer};
use core::hash::Hash;
use halo2_curves::group::ff::PrimeField;
use std::{collections::HashMap, iter::once, sync::Arc};

use crate::{backend::poly::multilinear::MultilinearPolynomial, Error};

/// Compute count polynomial, which is the number of occurrences of each table
/// entry `gx` in the input `fxs`.
/// Cost: linear in N.
pub(super) fn compute_count_in_input<F: PrimeField + Hash>(
    fxs: &[Arc<MultilinearPolynomial<F>>],
    gx: &Arc<MultilinearPolynomial<F>>,
) -> Result<Arc<MultilinearPolynomial<F>>, Error> {
    let start = start_timer!(|| "compute count in input");
    let g_evals = &gx.evals();

    let mut counter: HashMap<&F, usize> = g_evals
        .iter()
        .map(|g_eval| (g_eval, 0))
        .collect::<HashMap<&F, usize>>();
    for fx in fxs.into_iter() {
        let f_evals = &fx.evals();
        for f_eval in f_evals.iter() {
            if let Some(count) = counter.get_mut(f_eval) {
                *count += 1;
            } else {
                return Err(Error::InvalidSnark(
                    "failed to lookup fxs in gx".to_string(),
                ));
            }
        }
    }
    let count_evals = g_evals
        .iter()
        .map(|g_eval| match counter.remove(g_eval) {
            Some(count) => F::from(count as u64),
            None => F::zero(),
        })
        .collect();
    end_timer!(start);
    Ok(Arc::new(MultilinearPolynomial::new(count_evals)))
}

/// Compute the denominators and the entry fraction, where
/// denominators = [(fx + beta) | fx \in fxs] || [(gx + beta)]
/// h = \sum_{fx \in fxs} 1 / (fx + beta) - count / (gx + beta)
/// Cost: O(N)
pub(super) fn compute_h<F: PrimeField>(
    fxs: &[Arc<MultilinearPolynomial<F>>],
    gx: &Arc<MultilinearPolynomial<F>>,
    count: &Arc<MultilinearPolynomial<F>>,
    beta: &F,
) -> Result<Arc<MultilinearPolynomial<F>>, Error> {
    let start = start_timer!(|| "compute denominators and entry frac polynomial");
    let num_vars = gx.num_vars();
    // denominators = [(fx + beta) | fx \in fxs] || [(gx + beta)]
    let mut denominators = vec![];
    for poly in fxs.iter().chain(once(gx)) {
        if poly.num_vars() != num_vars {
            return Err(Error::InvalidSnark(
                "fxs have different number of variables".to_string(),
            ));
        }
        let mut denominator_evals = vec![];
        for eval in poly.iter() {
            denominator_evals.push(*eval + beta);
        }
        denominators.push(Arc::new(MultilinearPolynomial::new(denominator_evals)));
    }
    // compute h_evals = \sum_{f \in fs} 1/f - count / g
    let mut h_evals = vec![F::zero(); 1 << denominators[0].num_vars()];
    for denom in denominators.iter().take(denominators.len() - 1) {
        let denom_evals = &denom.evals();
        for (i, denom_eval) in denom_evals.iter().enumerate() {
            let frac = denom_eval.invert();
            if bool::from(frac.is_none()) {
                return Err(Error::InvalidSnark("(fx + beta) is zero".to_string()));
            }
            h_evals[i] += frac.unwrap();
        }
    }
    let last = denominators.last().unwrap();
    for (i, count_eval) in count.evals().iter().enumerate() {
        let frac = last.evals()[i].invert();
        if bool::from(frac.is_none()) {
            return Err(Error::InvalidSnark("(gx + beta) is zero".to_string()));
        }
        h_evals[i] -= *count_eval * frac.unwrap();
    }
    let h_poly = Arc::new(MultilinearPolynomial::new(h_evals));
    end_timer!(start);
    Ok(h_poly)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::poly::multilinear::{arc_mpoly, mpoly, MultilinearPolynomial};
    use halo2_curves::{bn256::Fr as F, group::ff::BatchInvert};
    use itertools::Itertools;
    use std::sync::Arc;

    #[test]
    fn test_compute_count_in_input() {
        let fxs = vec![arc_mpoly!(F, 1, 2, 3, 1), arc_mpoly!(F, 2, 2, 1, 3)];
        let gx = arc_mpoly!(F, 1, 2, 3, 4);
        let count = compute_count_in_input(&fxs, &gx).unwrap();
        assert_eq!(count.evals(), mpoly!(F, 3, 3, 2, 0).evals());
    }

    #[test]
    fn test_compute_denominators_and_h() {
        let fxs = vec![arc_mpoly!(F, 1, 2, 3, 1), arc_mpoly!(F, 2, 2, 1, 3)];
        let gx = arc_mpoly!(F, 1, 2, 3, 4);
        let count = compute_count_in_input(&fxs, &gx)
            .expect("Error occurred during compute_count_in_input");
        let beta = F::from(1u64);
        let h = compute_h(&fxs, &gx, &count, &beta).expect("Error occurred during compute_h");
        let denominators = vec![
            arc_mpoly!(F, 2, 3, 4, 2),
            arc_mpoly!(F, 3, 3, 2, 4),
            arc_mpoly!(F, 2, 3, 4, 5),
        ];
        let mut inv_denominators = denominators
            .iter()
            .map(|denom| {
                let mut inv_denom_evals = denom.evals().to_vec();
                inv_denom_evals.batch_invert();
                inv_denom_evals
            })
            .collect::<Vec<_>>();
        let last_term = inv_denominators.pop().unwrap();
        let last_term = last_term
            .iter()
            .zip(count.iter())
            .map(|(x, c)| -x * c)
            .collect_vec();
        let want_h = inv_denominators.iter().fold(last_term, |acc, x| {
            acc.iter().zip(x).map(|(a, b)| a + b).collect_vec()
        });

        assert_eq!(h.evals(), want_h);
    }
}
