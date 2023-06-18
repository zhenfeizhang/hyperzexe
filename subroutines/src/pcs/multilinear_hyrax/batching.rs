//! Sumcheck based batch opening and verify commitment.
// TODO: refactoring this code to somewhere else
// currently IOP depends on PCS because perm check requires commitment.
// The sumcheck based batch opening therefore cannot stay in the PCS repo --
// which creates a cyclic dependency.

use std::sync::Arc;

use crate::{
    pcs::{multilinear_hyrax::dense_mlpoly::EqPolynomial, PolynomialCommitmentScheme},
    HyraxCommitment, PCSError, HyraxProverParam, HyraxVerifierParam,
};
use arithmetic::{DenseMultilinearExtension};
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_std::{end_timer, start_timer, Zero};
use transcript::IOPTranscript;

/// A batch proof for polynomials open on a single point.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BatchProof<C, PCS>
where
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<C>,
{
    pub eval_groups: Vec<Vec<C::ScalarField>>,
    pub proof: PCS::Proof,
}

pub(crate) fn multi_open_single_point_internal<C, PCS>(
    prover_param: &PCS::ProverParam,
    polynomials: &[&[PCS::Polynomial]],
    point: &PCS::Point,
    eval_groups: &[&[C::ScalarField]],
    transcript: &mut IOPTranscript<C::ScalarField>,
) -> Result<BatchProof<C, PCS>, PCSError>
where
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<
        C,
        ProverParam = HyraxProverParam<C>,
        Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>,
        Point = Vec<C::ScalarField>,
    >,
{
    let open_timer = start_timer!(|| "multi open single point".to_string());
    // challenge extract point combination within a batch of polynomials
    let max_num_vars = prover_param.max_num_vars;

    // challenge to combine batches
    let eta = transcript.get_and_append_challenge(b"batch open RLC challenge")?;
    let mut poly_evals = vec![C::ScalarField::zero(); 1 << max_num_vars];
    for group in polynomials.iter() {
        let evs: Vec<C::ScalarField> = group.iter().flat_map(|p| p.evaluations.clone()).collect();
        for (i, ev) in evs.iter().enumerate() {
            poly_evals[i] = poly_evals[i] * eta + ev;
        }
    }
    // generate new point
    let mut new_point = point.clone();
    let point_high = transcript.get_and_append_challenge_vectors(
        b"batch open RLC challenge point",
        max_num_vars - point.len(),
    )?;
    new_point.extend(point_high.clone());
    let poly = Arc::new(DenseMultilinearExtension::from_evaluations_vec(
        max_num_vars,
        poly_evals,
    ));
    let (proof, eval) = PCS::open(prover_param, &poly, &new_point, transcript)?;

    // sanity check
    let rlc_coeff = EqPolynomial::new(point_high).evals();
    let mut expected_eval = C::ScalarField::zero();
    for group in eval_groups.iter() {
        expected_eval = expected_eval * eta;
        for (ev, co) in group.iter().zip(rlc_coeff.iter()) {
            expected_eval = expected_eval + *ev * *co;
        }
    }
    assert_eq!(eval, expected_eval);

    end_timer!(open_timer);

    Ok(BatchProof {
        eval_groups: eval_groups.iter().map(|g| g.to_vec()).collect(),
        proof,
    })
}

pub(crate) fn batch_verify_single_point_internal<C, PCS>(
    verifier_param: &PCS::VerifierParam,
    f_i_commitments: &[PCS::Commitment],
    point: &PCS::Point,
    batch_proof: &BatchProof<C, PCS>,
    transcript: &mut IOPTranscript<C::ScalarField>,
) -> Result<bool, PCSError>
where
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<
        C,
        VerifierParam = HyraxVerifierParam<C>,
        Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>,
        Point = Vec<C::ScalarField>,
        Commitment = HyraxCommitment<C>,
    >,
{
    let open_timer = start_timer!(|| "batch verification");
    let max_num_vars = verifier_param.max_num_vars;

    // randomly combine evaluations
    let eta = transcript.get_and_append_challenge(b"batch open RLC challenge")?;
    let point_high = transcript.get_and_append_challenge_vectors(
        b"batch open RLC challenge point",
        verifier_param.max_num_vars - point.len(),
    )?;
    let rlc_coeff = EqPolynomial::new(point_high.clone()).evals();
    let mut expected_eval = C::ScalarField::zero();
    for group in batch_proof.eval_groups.iter() {
        expected_eval = expected_eval * eta;
        for (eval, coeff) in group.iter().zip(rlc_coeff.iter()) {
            expected_eval = expected_eval + *eval * *coeff;
        }
    }

    // randomly combine commitments
    let combined_comm = {
        let (left_num_vars, _) = verifier_param.compute_factored_lens(max_num_vars);
        let mut combined_comm = vec![C::zero().into_projective(); 1 << left_num_vars];
        let eta_bigint = eta.into_repr();
        for f_i in f_i_commitments.iter() {
            for (i, c) in f_i.commitment.iter().enumerate() {
                combined_comm[i] = combined_comm[i].mul(eta_bigint) + (c.into_projective());
            }
        }
        HyraxCommitment {
            commitment: C::Projective::batch_normalization_into_affine(&combined_comm),
        }
    };

    // generate new point
    let mut new_point = point.clone();
    new_point.extend(point_high);
    let res = PCS::verify(
        verifier_param,
        &combined_comm,
        &new_point,
        &expected_eval,
        &batch_proof.proof,
        transcript,
    )?;

    end_timer!(open_timer);
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pcs::{
        multilinear_hyrax::{math::Math, nizk::DotProductProofGens},
        prelude::{MultilinearHyraxPCS},
    };
    use ark_bls12_381::{Fr, G1Affine as G1};
    use ark_ff::One;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::{rand::Rng, test_rng, vec::Vec, UniformRand};

    fn test_multi_open_helper<R: Rng>(
        srs: &DotProductProofGens<G1>,
        poly_groups: &[&[Arc<DenseMultilinearExtension<Fr>>]],
        rng: &mut R,
    ) -> Result<(), PCSError> {
        let num_vars = poly_groups[0][0].num_vars();
        let max_num_vars = poly_groups.iter().map(|g| g.len()).max().unwrap().next_power_of_two().log_2() + num_vars;
        let (ml_ck, ml_vk) = MultilinearHyraxPCS::<G1>::trim(srs, None, Some(max_num_vars))?;

        let point = (0..num_vars).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>();

        let evals: Vec<_> = poly_groups
            .iter()
            .map(|f| {
                f.iter()
                    .map(|f| f.evaluate(&point).unwrap())
                    .collect::<Vec<_>>()
            })
            .collect();
        let eval_slice = evals.iter().map(|e| e.as_slice()).collect::<Vec<_>>();

        let commitments: Vec<_> = poly_groups
            .iter()
            .map(|group| MultilinearHyraxPCS::multi_commit(&ml_ck, group).unwrap())
            .collect();

        let mut transcript = IOPTranscript::new("test transcript".as_ref());
        transcript.append_field_element("init".as_ref(), &Fr::zero())?;

        let mut batch_proof = multi_open_single_point_internal::<G1, MultilinearHyraxPCS<G1>>(
            &ml_ck,
            poly_groups,
            &point,
            eval_slice.as_slice(),
            &mut transcript,
        )?;

        // good path
        let mut transcript = IOPTranscript::new("test transcript".as_ref());
        transcript.append_field_element("init".as_ref(), &Fr::zero())?;
        assert!(batch_verify_single_point_internal(
            &ml_vk,
            commitments.as_slice(),
            &point,
            &batch_proof,
            &mut transcript
        )?);

        // bad path
        let mut wrong_evals = evals.clone();
        wrong_evals[0][0] = wrong_evals[0][0] + Fr::one();
        batch_proof.eval_groups = wrong_evals;

        let mut transcript = IOPTranscript::new("test transcript".as_ref());
        transcript.append_field_element("init".as_ref(), &Fr::zero())?;
        assert!(!batch_verify_single_point_internal(
            &ml_vk,
            commitments.as_slice(),
            &point,
            &batch_proof,
            &mut transcript
        )?);

        Ok(())
    }

    #[test]
    fn test_multi_open_internal() -> Result<(), PCSError> {
        let mut rng = test_rng();
        let srs = MultilinearHyraxPCS::gen_srs_for_testing(&mut rng, 20)?;

        for num_poly in 3..4 {
            for nv in 10..12 {
                let polys1: Vec<_> = (0..num_poly)
                    .map(|_| Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)))
                    .collect();
                let half = polys1.len() / 3;
                let poly_groups = polys1.split_at(half);
                
                test_multi_open_helper(&srs, &[poly_groups.0, poly_groups.1], &mut rng)?;
            }
        }

        Ok(())
    }
}
