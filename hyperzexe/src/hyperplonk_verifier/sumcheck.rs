use std::collections::HashMap;

use crate::{
    halo2_verifier::{loader::Loader, util::transcript::TranscriptRead},
    hyperplonk_verifier::sumcheck::sumcheck_round::ClassicSumcheckRoundProof,
};
use halo2_proofs::curves::CurveAffine;
use hyperplonk::backend::util::expression::{Expression, Query};

use self::sumcheck_round::SumcheckRoundVerifier;

use super::{util::poly::evaluate, Error};

pub mod classic;

pub trait SumcheckVerifier<C: CurveAffine, L: Loader<C>, SCR: SumcheckRoundVerifier<C, L>> {
    type Proof;
    type Output;

    fn read_proof<T>(num_vars: usize, degree: usize, transcript: &mut T) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        proof: &Self::Proof,
        expression: &Expression<L::LoadedScalar>,
        evals: &HashMap<Query, L::LoadedScalar>,
        challenges: &[L::LoadedScalar],
        ys: &[&[L::LoadedScalar]],
        sum: &L::LoadedScalar,
        num_vars: usize,
        degree: usize,
    ) -> Result<Self::Output, Error>;
}

pub trait SumcheckRoundVerifier<C: CurveAffine, L: Loader<C>> {
    type Proof: Clone + Debug;
    type Output: Clone + Debug;

    fn read_proof<T>(degree: usize, transcript: &mut T) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        proof: &Self::Proof,
        sum: &L::LoadedScalar,
        degree: usize,
        round: usize,
    ) -> Result<Self::Output, Error>;
}
