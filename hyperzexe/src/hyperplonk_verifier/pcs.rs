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
    type Proof: Clone + Debug;
    type Output: Clone + Debug;

    fn read_proof<T>(
        vk: &Self::VerifyingKey,
        queries: &[Query<C::Scalar>],
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        vk: &Self::VerifyingKey,
        commitments: &[Msm<C, L>],
        point: &L::LoadedScalar,
        queries: &[Query<C::Scalar, L::LoadedScalar>],
        proof: &Self::Proof,
    ) -> Result<Self::Output, Error>;
}
