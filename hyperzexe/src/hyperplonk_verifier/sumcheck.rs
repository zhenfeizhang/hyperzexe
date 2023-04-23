trait SumcheckVerifier<F: PrimeField, L: Loader<C>> {
    type Proof;
    type Output;

    fn read_proof<T>() -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify() -> Result<Self::Output, Error>;
}

struct SumcheckVerifier<F: PrimeField, L: Loader<C>> {}

impl<F: PrimeField, L: Loader<C>> SumcheckVerifier<F, L> {
    fn read_proof<T>() -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        unimplemented!()
    }

    fn verify() -> Result<Self::Output, Error> {
        unimplemented!()
    }
}
