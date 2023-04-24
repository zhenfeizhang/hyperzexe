mod sumcheck_round;

trait SumcheckVerifier<F: PrimeField, L: Loader<C>, SCR: SumcheckRoundVerifier<F, L>> {
    type VerifyingKey;
    type Proof;
    type Output;

    fn read_proof<T>(num_var: usize) -> Self::Proof
    where
        T: TranscriptRead<C, L>;

    fn verify(
        vk: &Self::VerifyingKey,
        expression: &Expression<L::LoadedScalar>,
        sum: &L::LoadedScalar,
        num_vars: usize,
        degree: usize,
    ) -> Result<Self::Output, Error>;
}

struct ClassicSumcheckVerifier<F: PrimeField, L: Loader<C>> {}

impl<F: PrimeField, L: Loader<C>> SumcheckVerifier<F, L> for ClassicSumcheckVerifier {
    type VerifyingKey = ();
    type Proof = SumcheckProof;
    type Output = ();
    fn read_proof<T>() -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        let round_proofs = Vec::new();
        for _ in 0..num_vars {
            round_proofs.push(SCR::read_proof::<T>());
        }
        ClassicSumcheckProof { round_proofs }
    }

    fn verify(
        _: &Self::VerifyingKey,
        expression: &Expression<L::LoadedScalar>,
        evals: &HashMap<Query, L::LoadedScalar>,
        challenges: &[L::LoadedScalar],
        ys: &[&[L::LoadedScalar]],
        x: &[L::LoadedScalar],
        sum: &L::LoadedScalar,
        num_vars: usize,
        degree: usize,
    ) -> Result<Self::Output, Error> {
        let mut round_sum = sum.clone();
        for i in 0..num_vars {
            round_sum = SCR::verify(_, round_proofs[i], &round_sum)?;
        }

        sum.loader().assert_eq(
            round_sum,
            evaluate(expression, num_vars, evals, challenges, ys, x),
        );
        Ok(())
    }
}

struct ClassicSumcheckProof<L> {
    round_proofs: Vec<&ClassicSumcheckRoundProof<L>>,
}
