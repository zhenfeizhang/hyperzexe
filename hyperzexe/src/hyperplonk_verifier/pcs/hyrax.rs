use std::marker::PhantomData;

use halo2_proofs::curves::CurveAffine;
use hyperplonk::backend::pcs::prelude::HyraxVerifierParam;

use crate::hyperplonk_verifier::pcs::MultiOpenScheme;

use super::OpenScheme;

mod bulletproof;
mod dot_product;

#[derive(Clone, Debug)]
pub struct Hyrax<C>(PhantomData<C>);

impl<C, L> OpenScheme<C, L> for Hyrax<C>
where
    C: CurveAffine,
    L: Loader<C>,
{
    type VerifyingKey = HyraxVerifierParam<C>;
    type Commitment = Vec<Msm<C, L>>;
    type Proof = HyraxOpenProof<C, L>;
    type Point = Vec<L::LoadedPoint>; // TODO: should be changed to a reference?
    type Output = ();

    fn read_proof<T>(
        vk: &Self::VerifyingKey,
        value: &L::LoadedScalar,
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        let dot_product_proof = DotProductProof::read_proof(vk, value, transcript);
        let blind_eval = value;

        Self::Proof {
            blind_eval,
            dot_product_proof,
        }
    }

    #[alow(non_snake_case)]
    fn verify(
        vk: &Self::VerifyingKey,
        commitments: &Self::Commitment,
        point: &Self::Point,
        proof: &Self::Proof,
    ) -> Result<Self::Output, Error> {
        // compute L and R
        let (_, right_num_vars) = gens.compute_factored_lens(point.len());
        let (L, R) = compute_factored_evals(point, right_num_vars);

        // compute a weighted sum of commitments and L
        let loader = value.loader();
        let C_LZ = L::multi_scalar_multiplication(L.iter().zip(comm.iter()));
        let C_Zr = L::multi_scalar_multiplication(
            iter::empty()
                .chain(iter::once(eval))
                .chain(iter::once(&self.blind_eval))
                .zip(gens.gens.gens_1.iter()),
        );
        self.proof
            .verify(R.len(), &gens.gens, transcript, &R, &C_LZ, &C_Zr)
    }
}

struct HyraxOpenProof<C, L> {
    blind_eval: L::LoadedScalar,
    dot_product_proof: DotProductProof<C, L>,
}

impl<C, L, OS> MultiOpenScheme<C, L, OS> for Hyrax<C>
where
    C: CurveAffine,
    L: Loader<C>,
    OS: OpenScheme<C, L>,
{
    type Proof = HyraxMultiOpenProof<C, L>;

    fn read_proof<T>(
        vk: &OS::VerifyingKey,
        num_var: usize,
        eval_groups: &[&L::LoadedScalar],
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        let eta = transcript.squeeze_challenge();
        let extra_challenges = transcript.squeeze_challenges(vk.max_num_vars - num_var);
        // randomly combine evaluations
        let rlc_coeff = eq_poly(&self.extra_challenges);
        let mut expected_eval = C::Scalar::zero();
        for group in batch_proof.eval_groups.iter() {
            expected_eval = expected_eval * eta;
            for (eval, coeff) in group.iter().zip(rlc_coeff.iter()) {
                expected_eval = expected_eval + *eval * *coeff;
            }
        }
        let proof = OS::read_proof(vk, expected_eval, transcript);
        Self::Proof {
            eval_groups,
            proof,
            eta,
            extra_challenges,
        }
    }

    fn verify(
        vk: &OS::VerifyingKey,
        commitments: &[&OS::Commitment],
        point: &L::LoadedScalar,
        batch_proof: &Self::Proof,
    ) -> Result<OS::Output, Error> {
        let max_num_vars = vk.max_num_vars;

        // randomly combine commitments
        let combined_comm = {
            let (left_num_vars, _) = verifier_param.compute_factored_lens(max_num_vars);
            let mut combined_comm: Vec<L::LoadedEcPoint> =
                vec![C::identity().into(); 1 << left_num_vars];
            for f_i in f_i_commitments.iter() {
                combined_comm.iter_mut().for_each(|c| *c = *c * eta);
                for (i, c) in f_i.iter().enumerate() {
                    combined_comm[i] = combined_comm[i] + c;
                }
            }
            combined_comm
        };

        // generate new point
        let mut new_point = point.clone();
        new_point.extend(point_high);
        let res = OS::verify(
            verifier_param,
            &combined_comm,
            &new_point,
            &batch_proof.proof,
        )?;
    }
}
struct HyraxMultiOpenProof<C, L> {
    pub eval_groups: Vec<Vec<L::LoadedScalar>>,
    pub proof: HyraxOpenProof<C, L>,
    // random linear combination of commitments
    pub eta: L::LoadedScalar,
    // extra challenges to construct a new point
    pub extra_challenges: Vec<L::LoadedScalar>,
}
