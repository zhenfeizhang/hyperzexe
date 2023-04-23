use std::marker::PhantomData;

use halo2_proofs::curves::CurveAffine;
use hyperplonk::backend::pcs::prelude::HyraxVerifierParam;

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
    type Proof = ();
    type Output = ();

    fn read_proof<T>(
        _vk: &Self::VerifyingKey,
        _queries: &[Query<C::Scalar>],
        _transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        unimplemented!()
    }

    fn verify(
        _vk: &Self::VerifyingKey,
        _commitments: &[Msm<C, L>],
        _point: &L::LoadedScalar,
        _queries: &[Query<C::Scalar, L::LoadedScalar>],
        _proof: &Self::Proof,
    ) -> Result<Self::Output, Error> {
        unimplemented!()
    }
}

#[allow(non_snake_case)]
struct HyraxBatchProof<C, L> {
    pub eval_groups: Vec<Vec<L::LoadedScalar>>,
    // bulletproof
    dot_product_proof: DotProductProof<C, L>,
}
