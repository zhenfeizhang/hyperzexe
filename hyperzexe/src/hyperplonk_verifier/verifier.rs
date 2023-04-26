use halo2_proofs::curves::CurveAffine;
use std::fmt::Debug;

use crate::{
    halo2_verifier::{loader::Loader, util::transcript::TranscriptRead},
    hyperplonk_verifier::Protocol,
};

use super::{
    pcs::{MultiOpenScheme, OpenScheme},
    sumcheck::{sumcheck_round::SumcheckRoundVerifier, SumcheckVerifier},
    Error,
};

mod hyperplonk;

pub trait HyperPlonkVerifier<C, L, SRC, SC, OS, MOS>
where
    C: CurveAffine,
    L: Loader<C>,
    SRC: SumcheckRoundVerifier<C, L>,
    SC: SumcheckVerifier<C, L, SRC>,
    OS: OpenScheme<C, L>,
    MOS: MultiOpenScheme<C, L, OS>,
{
    type VerifyingKey: Clone + Debug;
    type Proof: Clone + Debug;
    type Output: Clone + Debug;

    fn read_proof<T>(
        vk: &Self::VerifyingKey,
        protocol: &Protocol<C, L>,
        instances: &[Vec<L::LoadedScalar>],
        transcript: &mut T,
    ) -> Result<Self::Proof, Error>
    where
        T: TranscriptRead<C, L>;

    fn verify(
        vk: &Self::VerifyingKey,
        protocol: &Protocol<C, L>,
        instances: &[Vec<L::LoadedScalar>],
        proof: &Self::Proof,
    ) -> Result<Self::Output, Error>;
}
