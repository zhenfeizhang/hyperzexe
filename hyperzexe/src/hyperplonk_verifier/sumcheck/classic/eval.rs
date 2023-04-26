use std::{iter, marker::PhantomData};

use halo2_proofs::curves::CurveAffine;

use crate::{
    halo2_verifier::{loader::Loader, util::transcript::TranscriptRead},
    hyperplonk_verifier::{
        util::poly::{barycentric_interpolate, barycentric_weights},
        Error,
    },
};
use std::fmt::Debug;

pub struct ClassicSumcheckRoundVerifier<C: CurveAffine, L: Loader<C>> {
    _marker: PhantomData<(C, L)>,
}

impl<C: CurveAffine, L: Loader<C>> SumcheckRoundVerifier<C, L>
    for ClassicSumcheckRoundVerifier<C, L>
{
    type Proof = ClassicSumcheckRoundProof<C, L>;
    type Output = L::LoadedScalar;
    fn read_proof<T>(degree: usize, transcript: &mut T) -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        let msg = transcript.read_n_scalars(degree + 1)?;
        let challenge = transcript.squeeze_challenge();
        ClassicSumcheckRoundProof { challenge, msg }
    }

    fn verify(
        proof: &Self::Proof,
        sum: &L::LoadedScalar,
        degree: usize,
        round: usize,
    ) -> Result<Self::Output, Error> {
        let loader = sum.loader();
        loader.assert_eq(&proof.msg.sum(), sum);
        // TODO: precompute and store it.
        let native_zero = C::Scalar::zero();
        let native_one = C::Scalar::one();
        let points = iter::successors(Some(native_zero), move |state| Some(native_one + state))
            .take(degree + 1)
            .map(|x| loader.load_const(x))
            .collect::<Vec<_>>();
        let weight = barycentric_weights(&points);
        let sum = barycentric_interpolate(&weight, &points, &proof.msg, &proof.challenge);
        Ok(sum)
    }
}

pub struct ClassicSumcheckRoundProof<C: CurveAffine, L: Loader<C>> {
    pub challenge: L::LoadedScalar,
    pub msg: Vec<L::LoadedScalar>,
}
