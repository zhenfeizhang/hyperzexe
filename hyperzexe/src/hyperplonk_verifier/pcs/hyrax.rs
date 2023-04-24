pub trait MultiOpenScheme<C, L>: Hyrax<C, L>
where
    C: CurveAffine,
    L: Loader<C>,
{
    type SuccinctVerifyingKey: Clone + Debug;
    type Proof: Clone + Debug;

    fn read_proof<T>(
        svk: &Self::SuccinctVerifyingKey,
        queries: &[Query<C::Scalar>],
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn succinct_verify(
        svk: &Self::SuccinctVerifyingKey,
        commitments: &[Msm<C, L>],
        point: &L::LoadedScalar,
        queries: &[Query<C::Scalar, L::LoadedScalar>],
        proof: &Self::Proof,
    ) -> Self::Accumulator;
}
