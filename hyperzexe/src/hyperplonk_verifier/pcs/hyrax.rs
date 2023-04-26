use std::{iter, marker::PhantomData};

use halo2_proofs::curves::CurveAffine;
use hyperplonk::backend::pcs::prelude::HyraxVerifierParam;

use crate::{
    halo2_verifier::{
        loader::Loader,
        util::{msm::Msm, transcript::TranscriptRead},
    },
    hyperplonk_verifier::{
        pcs::MultiOpenScheme,
        util::poly::{compute_factored_evals, eq_eval},
        Error,
    },
};

use self::dot_product::DotProductProof;

use super::OpenScheme;

mod bulletproof;
mod dot_product;

#[derive(Clone, Debug)]
pub struct Hyrax<C>(PhantomData<C>);

impl<'a, C, L> OpenScheme<C, L> for Hyrax<C>
where
    C: CurveAffine,
    L: Loader<C>,
{
    type VerifyingKey = HyraxVerifierParam<C>;
    type Commitment = Vec<L::LoadedEcPoint>;
    type Proof = HyraxOpenProof<C, L>;
    type Point = Vec<L::LoadedScalar>; // TODO: should be changed to a reference?
    type Output = ();

    fn read_proof<T>(blind: &L::LoadedScalar, transcript: &mut T, n: usize) -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        let dot_product_proof = DotProductProof::read_proof(n, transcript);
        let blind_eval = blind;

        Self::Proof {
            blind_eval,
            dot_product_proof,
        }
    }

    #[allow(non_snake_case)]
    fn verify(
        vk: &Self::VerifyingKey,
        commitments: &Self::Commitment,
        point: &Self::Point,
        value: &L::LoadedScalar,
        proof: &Self::Proof,
    ) -> Result<Self::Output, Error> {
        // compute L and R
        let (_, right_num_vars) = vk.compute_factored_lens(point.len());
        let (L, R) = compute_factored_evals(point, right_num_vars);

        // compute a weighted sum of commitments and L
        let loader = proof.blind_eval.loader();
        let C_LZ = L::multi_scalar_multiplication(L.iter().zip(commitments.iter()));
        let C_Zr = L::multi_scalar_multiplication(
            iter::empty()
                .chain(iter::once(value))
                .chain(iter::once(&proof.blind_eval))
                .zip([vk.gens.gens_1.G, vk.gens.gens_1.h].iter()),
        );
        proof.dot_product_proof.verify(&vk.gens, &R, &C_LZ, &C_Zr)
    }
}

struct HyraxOpenProof<C: CurveAffine, L: Loader<C>> {
    blind_eval: L::LoadedScalar,
    dot_product_proof: DotProductProof<C, L>,
}

impl<C, L, OS> MultiOpenScheme<C, L, OS> for Hyrax<C>
where
    C: CurveAffine,
    L: Loader<C>,
    OS: OpenScheme<C, L, VerifyingKey = HyraxVerifierParam<C>>,
{
    type VerifyingKey = OS::VerifyingKey;
    type Commitment = OS::Commitment;
    type Point = OS::Point;
    type Proof = HyraxMultiOpenProof<C, L>;
    type Output = OS::Output;

    fn read_proof<T>(
        vk: &Self::VerifyingKey,
        num_var: usize,
        eval_groups: &[&L::LoadedScalar],
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        let eta = transcript.squeeze_challenge();
        let extra_challenges = transcript.squeeze_challenges(vk.max_num_vars - num_var);
        let blind = eval_groups[0].loader().load_zero();
        let proof = OS::read_proof(vk, blind, transcript);
        Self::Proof {
            eval_groups,
            proof,
            eta,
            extra_challenges,
        }
    }

    fn verify(
        vk: &Self::VerifyingKey,
        commitments: &[&Self::Commitment],
        point: &Self::Point,
        batch_proof: &Self::Proof,
    ) -> Result<Self::Output, Error> {
        let max_num_vars = vk.max_num_vars;

        // randomly combine commitments
        let combined_comm = {
            let (left_num_vars, _) = vk.compute_factored_lens(max_num_vars);
            let mut combined_comm: Vec<L::LoadedEcPoint> =
                vec![C::identity().into(); 1 << left_num_vars];
            for commitment in commitments.iter() {
                combined_comm
                    .iter_mut()
                    .for_each(|c| *c = *c * batch_proof.eta);
                for (i, c) in commitment.iter().enumerate() {
                    combined_comm[i] = combined_comm[i] + c;
                }
            }
            combined_comm
        };

        // randomly combine evaluations
        let rlc_coeff = eq_eval(&batch_proof.extra_challenges);
        let loader = batch_proof.eval_groups[0].loader();
        let mut expected_eval = loader.load_const(C::Scalar::zero());
        for group in batch_proof.eval_groups.iter() {
            expected_eval = expected_eval * batch_proof.eta;
            for (eval, coeff) in group.iter().zip(rlc_coeff.iter()) {
                expected_eval = expected_eval + *eval * *coeff;
            }
        }

        // generate new point
        let mut new_point = point.clone();
        new_point.extend(batch_proof.extra_challenges);
        let res = OS::verify(
            vk,
            &combined_comm,
            &new_point,
            &expected_eval,
            &batch_proof.proof,
        )?;
    }

    fn compute_comm_len(vk: &Self::VerifyingKey, total_len: usize) -> usize {
        total_len >> vk.right_num_vars
    }
}
struct HyraxMultiOpenProof<C: CurveAffine, L: Loader<C>> {
    pub eval_groups: Vec<Vec<L::LoadedScalar>>,
    pub proof: HyraxOpenProof<C, L>,
    // random linear combination of commitments
    pub eta: L::LoadedScalar,
    // extra challenges to construct a new point
    pub extra_challenges: Vec<L::LoadedScalar>,
}
