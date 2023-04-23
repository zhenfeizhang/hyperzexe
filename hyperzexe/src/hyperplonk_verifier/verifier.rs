use halo2_proofs::curves::CurveAffine;
use std::fmt::Debug;

use crate::{hyperplonk_verifier::Protocol, halo2_verifier::util::transcript::TranscriptRead};

mod hyperplonk;

pub trait HyperPlonkVerifier<C, L, PCS>
where
    C: CurveAffine,
    L: Loader<C>,
{
    type Proof: Clone + Debug;

    fn read_proof<T>(
        svk: &PCS::SuccinctVerifyingKey,
        protocol: &Protocol<C, L>,
        instances: &[Vec<L::LoadedScalar>],
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        svk: &PCS::SuccinctVerifyingKey,
        protocol: &Protocol<C, L>,
        instances: &[Vec<L::LoadedScalar>],
        proof: &Self::Proof,
    ) -> PCS::Output;
}