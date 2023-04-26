use halo2_proofs::{curves::CurveAffine, ff::PrimeField};
use std::fmt::Debug;

use crate::halo2_verifier::{loader::Loader, util::transcript::TranscriptRead};

use super::Error;

pub mod hyrax;

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
        blind: &L::LoadedScalar,
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        vk: &Self::VerifyingKey,
        commitments: &Self::Commitment,
        point: &Self::Point,
        value: &L::LoadedScalar,
        proof: &Self::Proof,
    ) -> Result<Self::Output, Error>;
}

pub trait MultiOpenScheme<C, L, OS>
where
    C: CurveAffine,
    L: Loader<C>,
    OS: OpenScheme<C, L>,
{
    type VerifyingKey: Clone + Debug;
    type Commitment: Clone + Debug;
    type Point: Clone + Debug;
    type Proof: Clone + Debug;
    type Output: Clone + Debug;

    fn read_proof<T>(
        vk: &Self::VerifyingKey,
        num_vars: usize,
        eval_groups: &[&L::LoadedScalar],
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        vk: &Self::VerifyingKey,
        commitments: &[&Self::Commitment],
        point: &Self::Point,
        batch_proof: &Self::Proof,
    ) -> Result<Self::Output, Error>;

    fn compute_comm_len(vk: &Self::VerifyingKey, total_len: usize) -> usize;
}
