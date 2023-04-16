//! Sumcheck based batch opening and verify commitment.
// TODO: refactoring this code to somewhere else
// currently IOP depends on PCS because perm check requires commitment.
// The sumcheck based batch opening therefore cannot stay in the PCS repo --
// which creates a cyclic dependency.

use std::sync::Arc;

use halo2_curves::{
    group::{ff::Field, Curve},
    CurveAffine,
};

use crate::{
    backend::{
        pcs::PolynomialCommitmentScheme,
        poly::multilinear::MultilinearPolynomial,
        util::{end_timer, start_timer, transcript::TranscriptWrite},
    },
    Error,
};

use super::{HyraxProverParam, HyraxVerifierParam};

/// A batch proof for polynomials open on a single point.
#[derive(Clone, Debug)]
pub struct HyraxBatchProof<F, PCS>
where
    F: Field,
    PCS: PolynomialCommitmentScheme<F>,
{
    pub eval_groups: Vec<Vec<F>>,
    pub proof: PCS::Proof,
}

pub(crate) fn multi_open_single_point_internal<C, PCS>(
    prover_param: &PCS::ProverParam,
    polynomials: &[&[PCS::Polynomial]],
    point: &PCS::Point,
    eval_groups: &[&[C::Scalar]],
    transcript: &mut impl TranscriptWrite<C, C::Scalar>,
) -> Result<HyraxBatchProof<C::Scalar, PCS>, Error>
where
    C: CurveAffine,
    PCS: PolynomialCommitmentScheme<
        C::Scalar,
        Curve = C,
        ProverParam = HyraxProverParam<C>,
        Polynomial = Arc<MultilinearPolynomial<C::Scalar>>,
        Point = Vec<C::Scalar>,
    >,
{
    let open_timer = start_timer(|| "multi open single point".to_string());
    // challenge extract point combination within a batch of polynomials
    let max_num_vars = prover_param.max_num_vars;

    // challenge to combine batches
    let eta = transcript.squeeze_challenge();
    let mut poly_evals = vec![C::Scalar::zero(); 1 << max_num_vars];
    for group in polynomials.iter() {
        poly_evals.iter_mut().for_each(|e| *e = *e * eta);
        let evs: Vec<C::Scalar> = group.iter().flat_map(|p| p.evals().to_vec()).collect();
        for (i, ev) in evs.iter().enumerate() {
            poly_evals[i] = poly_evals[i] + ev;
        }
    }
    // generate new point
    let mut new_point = point.clone();
    let point_high = transcript.squeeze_challenges(max_num_vars - point.len());
    new_point.extend(point_high.clone());
    let poly = Arc::new(MultilinearPolynomial::new(poly_evals));
    let (proof, eval) = PCS::open(prover_param, &poly, &new_point, transcript)?;

    // sanity check
    let rlc_coeff = MultilinearPolynomial::eq_xy(&point_high).into_evals();
    let mut expected_eval = C::Scalar::zero();
    for group in eval_groups.iter() {
        expected_eval = expected_eval * eta;
        for (ev, co) in group.iter().zip(rlc_coeff.iter()) {
            expected_eval = expected_eval + *ev * *co;
        }
    }
    assert_eq!(eval, expected_eval);

    end_timer(open_timer);

    Ok(HyraxBatchProof {
        eval_groups: eval_groups.iter().map(|g| g.to_vec()).collect(),
        proof,
    })
}

pub(crate) fn batch_verify_single_point_internal<C, PCS>(
    verifier_param: &PCS::VerifierParam,
    f_i_commitments: &[&PCS::Commitment],
    point: &PCS::Point,
    batch_proof: &HyraxBatchProof<C::Scalar, PCS>,
    transcript: &mut impl TranscriptWrite<C, C::Scalar>,
) -> Result<bool, Error>
where
    C: CurveAffine,
    PCS: PolynomialCommitmentScheme<
        C::Scalar,
        Curve = C,
        VerifierParam = HyraxVerifierParam<C>,
        Polynomial = Arc<MultilinearPolynomial<C::Scalar>>,
        Point = Vec<C::Scalar>,
        Commitment = Vec<C>,
    >,
{
    let open_timer = start_timer(|| "batch verification");
    let max_num_vars = verifier_param.max_num_vars;

    // randomly combine evaluations
    let eta = transcript.squeeze_challenge();
    let point_high = transcript.squeeze_challenges(verifier_param.max_num_vars - point.len());
    let rlc_coeff = MultilinearPolynomial::eq_xy(&point_high).into_evals();
    let mut expected_eval = C::Scalar::zero();
    for group in batch_proof.eval_groups.iter() {
        expected_eval = expected_eval * eta;
        for (eval, coeff) in group.iter().zip(rlc_coeff.iter()) {
            expected_eval = expected_eval + *eval * *coeff;
        }
    }

    // randomly combine commitments
    let combined_comm = {
        let (left_num_vars, _) = verifier_param.compute_factored_lens(max_num_vars);
        let mut combined_comm: Vec<C::CurveExt> = vec![C::identity().into(); 1 << left_num_vars];
        for f_i in f_i_commitments.iter() {
            combined_comm.iter_mut().for_each(|c| *c = *c * eta);
            for (i, c) in f_i.iter().enumerate() {
                combined_comm[i] = combined_comm[i] + c;
            }
        }
        let mut commitment = vec![C::identity(); combined_comm.len()];
        C::Curve::batch_normalize(&combined_comm, &mut commitment);
        commitment
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

    end_timer(open_timer);
    Ok(res)
}

#[cfg(test)]
mod tests {
    use halo2_curves::bn256::{Fr as F, G1Affine as C};
    use itertools::Itertools;

    use super::*;
    use crate::backend::{
        pcs::{multilinear_hyrax::nizk::DotProductProofGens, prelude::MultilinearHyraxPCS},
        util::{
            test::std_rng,
            transcript::Keccak256Transcript,
        }, poly::multilinear::arc_mpoly,
    };

    fn test_multi_open_helper(
        srs: &DotProductProofGens<C>,
        poly_groups: &[&[Arc<MultilinearPolynomial<F>>]],
        point: &Vec<F>,
        evals: &[&[F]],
    ) -> Result<(), Error> {
        let num_vars = poly_groups[0][0].num_vars();
        let max_num_batches = poly_groups.iter().map(|g| g.len()).max().unwrap();
        let (ml_ck, ml_vk) =
            MultilinearHyraxPCS::<C>::trim(srs, None, Some(num_vars), max_num_batches)?;

        let commitments: Vec<_> = poly_groups
            .iter()
            .map(|group| MultilinearHyraxPCS::multi_commit(&ml_ck, *group).unwrap())
            .collect();

        let mut transcript = Keccak256Transcript::<Vec<u8>>::default();
        let mut batch_proof = multi_open_single_point_internal::<C, MultilinearHyraxPCS<C>>(
            &ml_ck,
            poly_groups,
            point,
            evals,
            &mut transcript,
        )?;

        // good path
        let commitment_slice = commitments.iter().collect_vec();
        let mut transcript = Keccak256Transcript::<Vec<u8>>::default();
        assert!(batch_verify_single_point_internal(
            &ml_vk,
            commitment_slice.as_slice(),
            point,
            &batch_proof,
            &mut transcript
        )?);

        // bad path
        let mut wrong_evals = evals.iter().map(|eval| eval.to_vec()).collect_vec();
        wrong_evals[0][0] = wrong_evals[0][0] + F::one();
        batch_proof.eval_groups = wrong_evals;

        let mut transcript = Keccak256Transcript::<Vec<u8>>::default();
        assert!(!batch_verify_single_point_internal(
            &ml_vk,
            commitment_slice.as_slice(),
            point,
            &batch_proof,
            &mut transcript
        )?);

        Ok(())
    }

    #[test]
    fn test_multi_open_internal() -> Result<(), Error> {
        let mut rng = std_rng();
        let srs = MultilinearHyraxPCS::<C>::setup(&mut rng, 20)?;

        let poly_group1 = vec![
            arc_mpoly!(F, 2, 8, 3, 5, 4, 7, 9, 10),
            arc_mpoly!(F, 1, 2, 2, 3, 7, 9, 1, 9),
        ];
        let poly_group2 = vec![
            arc_mpoly!(F, 1, 4, 1, 5, 5, 7, 3, 10),
        ];

        let evals = vec![vec![F::from(12078u64), F::from(31649u64)], vec![F::from(22138u64)]];
        let eval_slices = evals.iter().map(|e| e.as_slice()).collect::<Vec<_>>();
        let point = vec![F::from(11u64), F::from(17u64), F::from(31u64)];
        test_multi_open_helper(&srs, &[&poly_group1, &poly_group2], &point, &eval_slices)?;

        Ok(())
    }
}
