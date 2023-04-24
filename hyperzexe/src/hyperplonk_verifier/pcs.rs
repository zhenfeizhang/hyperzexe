use halo2_proofs::{curves::CurveAffine, ff::PrimeField};
use std::fmt::Debug;

use crate::halo2_verifier::util::transcript::TranscriptRead;

pub mod hyrax;

#[derive(Clone, Debug)]
pub struct Query<F: PrimeField, T = ()> {
    pub poly: usize,
    pub shift: F,
    pub eval: T,
}

pub trait OpenScheme<C, L>
where
    C: CurveAffine,
    L: Loader<C>,
{
    type VerifyingKey: Clone + Debug;
    type Commitment: Clone + Debug;
    type Proof: Clone + Debug;
    type Point: Clone + Debug;
    type Output: Clone + Debug;

    fn read_proof<T>(
        vk: &Self::VerifyingKey,
        value: &L::LoadedScalar,
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        vk: &Self::VerifyingKey,
        commitments: &Self::Commitment,
        point: &Self::Point,
        proof: &Self::Proof,
    ) -> Result<Self::Output, Error>;
}

pub trait MultiOpenScheme<C, L, OS>
where
    C: CurveAffine,
    L: Loader<C>,
    OS: OpenScheme<C, L>,
{
    type Proof: Clone + Debug;
    fn read_proof<T>(
        vk: &OS::VerifyingKey,
        num_vars: usize,
        eval_groups: &[&L::LoadedScalar],
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        vk: &OS::VerifyingKey,
        commitments: &[&OS::Commitment],
        point: &L::LoadedScalar,
        batch_proof: &Self::Proof,
    ) -> Result<OS::Output, Error>;
}
