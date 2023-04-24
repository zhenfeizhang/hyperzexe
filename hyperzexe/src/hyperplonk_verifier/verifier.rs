use ::hyperplonk::halo2_curves::pasta::pallas::Scalar;
use halo2_proofs::curves::CurveAffine;
use std::fmt::Debug;

use crate::{halo2_verifier::util::transcript::TranscriptRead, hyperplonk_verifier::Protocol};

use super::ScalarPoint;

mod hyperplonk;

pub trait HyperPlonkVerifier<C, PCS>
where
    C: CurveAffine,
{
    type Proof: Clone + Debug;

    fn read_proof<T, L>(
        vk: &PCS::VerifyingKey,
        protocol: &Protocol<C>,
        instances: &[Vec<ScalarPoint<C>>],
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        vk: &PCS::VerifyingKey,
        protocol: &Protocol<C>,
        instances: &[Vec<ScalarPoint<C>>],
        proof: &Self::Proof,
    ) -> PCS::Output;
}
