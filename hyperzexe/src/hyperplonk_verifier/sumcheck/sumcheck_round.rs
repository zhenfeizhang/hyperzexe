trait SumcheckRoundVerifier<F: PrimeField, L: Loader<C>> {
    type VerifyingKey;
    type Proof;
    type Output;

    fn read_proof<T>() -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify() -> Result<Self::Output, Error>;
}

struct ClassicSumcheckRoundVerifier<F: PrimeField, L: Loader<C>> {}

impl<F: PrimeField, L: Loader<C>> ClassicSumcheckRoundVerifier<F, L> {
    type VerifyingKey = ();
    type Proof = ClassicSumcheckRoundProof;
    type Output = L::LoadedScalar;
    fn read_proof<T>(
        degree: usize,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        let msg = transcript.read_n_scalars(degree + 1)?;
        challenge = transcript.squeeze_challenge();
        ClassicSumcheckRoundProof { challenge, msg }
    }

    fn verify(
        _: &Self::VerifyingKey,
        proof: &Self::Proof,
        sum: &L::LoadedScalar,
    ) -> Result<Self::Output, Error> {
        if sum != msg.sum() {
            let msg = if round == 0 {
                format!("Expect sum {sum:?} but get {:?}", msg.sum())
            } else {
                format!("Consistency failure at round {round}")
            };
            return Err(Error::InvalidSumcheck(msg));
        }
        let loader = sum.loader();
        let zero = loader.load_zero();
        let one = loader.load_one();
        // TODO: precompute and store it.
        let points = iter::successors(Some(zero), move |state| Some(one + state))
            .take(degree + 1)
            .collect_vec();
        let weight = barycentric_weights(points);
        let sum = barycentric_interpolate<L>(&weight, &points, &msg, &challenge);
        Ok(sum)
    }
}

struct ClassicSumcheckRoundProof<L> {
    challenge: L::LoadedScalar,
    msg: Vec<L::LoadedScalar>,
}
