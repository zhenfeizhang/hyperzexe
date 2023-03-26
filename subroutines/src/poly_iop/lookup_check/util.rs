use crate::poly_iop::{errors::PolyIOPErrors, structs::IOPProof, zero_check::ZeroCheck, PolyIOP};
use arithmetic::VirtualPolynomial;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use std::{collections::HashMap, iter::once, sync::Arc};
use transcript::IOPTranscript;

/// Compute the number of occurrences of each table entry `gx` in the input
/// `fxs`. Cost: linear in N.
pub(super) fn compute_count_in_input<F: PrimeField>(
    fxs: &[Arc<DenseMultilinearExtension<F>>],
    gx: &Arc<DenseMultilinearExtension<F>>,
) -> Result<Arc<DenseMultilinearExtension<F>>, PolyIOPErrors> {
    let start = start_timer!(|| "compute count in input");
    let g_evals = &gx.evaluations;
    let num_vars = gx.num_vars;

    let mut counter: HashMap<&F, usize> = g_evals
        .iter()
        .map(|g_eval| (g_eval, 0))
        .collect::<HashMap<&F, usize>>();
    for fx in fxs.iter() {
        let f_evals = &fx.evaluations;
        for f_eval in f_evals.iter() {
            if let Some(count) = counter.get_mut(f_eval) {
                *count += 1;
            } else {
                return Err(PolyIOPErrors::InvalidParameters(
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
    Ok(Arc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars,
        count_evals,
    )))
}

/// Compute the denominators and the entry fraction, where
/// denominators = [(fx + beta) | fx \in fxs] || [(gx + beta)]
/// h = \sum_{fx \in fxs} 1 / (fx + beta) - count / (gx + beta)
/// Cost: O(N)
pub(super) fn compute_denominators_and_entry_frac<F: PrimeField>(
    fxs: &[Arc<DenseMultilinearExtension<F>>],
    gx: &Arc<DenseMultilinearExtension<F>>,
    count: &Arc<DenseMultilinearExtension<F>>,
    beta: F,
) -> Result<
    (
        Vec<Arc<DenseMultilinearExtension<F>>>,
        Arc<DenseMultilinearExtension<F>>,
    ),
    PolyIOPErrors,
> {
    let start = start_timer!(|| "compute denominators and entry frac polynomial");
    let num_vars = gx.num_vars;
    // denominators = [(fx + beta) | fx \in fxs] || [(gx + beta)]
    let mut denominators = vec![];
    for poly in fxs.iter().chain(once(gx)) {
        if poly.num_vars != num_vars {
            return Err(PolyIOPErrors::InvalidParameters(
                "fxs have different number of variables".to_string(),
            ));
        }
        let mut denominator_evals = vec![];
        for eval in poly.iter() {
            denominator_evals.push(*eval + beta);
        }
        denominators.push(Arc::new(DenseMultilinearExtension::from_evaluations_vec(
            num_vars,
            denominator_evals,
        )));
    }
    // compute h_evals = \sum_{f \in fs} 1/f - count / g
    let mut h_evals = vec![F::zero(); denominators[0].evaluations.len()];
    for denom in denominators.iter().take(denominators.len() - 1) {
        let denom_evals = &denom.evaluations;
        for (i, denom_eval) in denom_evals.iter().enumerate() {
            let frac = denom_eval.inverse();
            if frac.is_none() {
                return Err(PolyIOPErrors::InvalidParameters(
                    "(fx + beta) is zero".to_string(),
                ));
            }
            h_evals[i] += frac.unwrap();
        }
    }
    let last = denominators.last().unwrap();
    for (i, count_eval) in count.evaluations.iter().enumerate() {
        let frac = last.evaluations[i].inverse();
        if frac.is_none() {
            return Err(PolyIOPErrors::InvalidParameters(
                "(gx + beta) is zero".to_string(),
            ));
        }
        h_evals[i] -= *count_eval * frac.unwrap();
    }
    let h_poly = Arc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars, h_evals,
    ));
    end_timer!(start);
    Ok((denominators, h_poly))
}

/// Generate the zerocheck proof for the virtual polynomial
///    (h(x)(gx(x) + beta) + count(x))\prod_{fx \in fxs}(fx(x) + beta)
///    - gx(x)\sum_(fx \in fxs) \prod_{fx' \neq fx}(fx'(x) + beta)
/// Returns proof.
/// - input: denominators = [f1(x) + beta, .., fk(x) + beta, g(x) + beta], h =
///   \sum_{fx \in fxs} 1 / (fx + beta) - count / (gx + beta)
///
/// Cost: O(N)
pub(super) fn prove_zero_check<F: PrimeField>(
    denominators: &[Arc<DenseMultilinearExtension<F>>],
    h: &Arc<DenseMultilinearExtension<F>>,
    count: &Arc<DenseMultilinearExtension<F>>,
    transcript: &mut IOPTranscript<F>,
) -> Result<(IOPProof<F>, VirtualPolynomial<F>), PolyIOPErrors> {
    let start = start_timer!(|| "zerocheck in lookup check");
    let num_vars = denominators[0].num_vars;

    let mut list_l = denominators.to_vec();
    let mut list_r: Vec<Arc<DenseMultilinearExtension<F>>> = Vec::new();

    let mut q_x = VirtualPolynomial::new(num_vars);
    // q(x) = (\prod_{fx \in fxs} fx(x)) * gx(x) * h(x)
    q_x.add_mle_list([list_l.clone(), vec![h.clone()]].concat(), F::one())?;

    // + \prod_{fx \in fxs}(fx(x) + beta) * count(x)
    let curr = list_l.pop().unwrap();
    q_x.add_mle_list([list_l.clone(), vec![count.clone()]].concat(), F::one())?;
    list_r.push(curr);

    // - gx(x)\sum_{fx' \in fxs}\prod_{fx \in fxs}(fx(x) + beta) / (fx'(x) + beta)
    for _ in 1..denominators.len() {
        let curr = list_l.pop().unwrap();
        q_x.add_mle_list([list_l.clone(), list_r.clone()].concat(), -F::one())?;
        list_r.push(curr);
    }

    let iop_proof = <PolyIOP<F> as ZeroCheck<F>>::prove(&q_x, transcript)?;
    end_timer!(start);
    Ok((iop_proof, q_x))
}
