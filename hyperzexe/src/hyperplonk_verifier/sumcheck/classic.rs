mod eval;

pub struct ClassicSumcheckVerifier<C: CurveAffine, L: Loader<C>> {
    _marker: std::marker::PhantomData<(C, L)>,
}

impl<C: CurveAffine, L: Loader<C>, SCR: SumcheckRoundVerifier<C, L>> SumcheckVerifier<C, L, SCR>
    for ClassicSumcheckVerifier<C, L>
{
    type Proof = ClassicSumcheckProof<C, L>;
    type Output = Vec<L::LoadedScalar>;
    fn read_proof<T>(num_vars: usize, degree: usize, transcript: &mut T) -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        let round_proofs = Vec::new();
        for _ in 0..num_vars {
            round_proofs.push(SCR::read_proof::<T>(degree, transcript));
        }
        ClassicSumcheckProof { round_proofs }
    }

    fn verify(
        proof: &Self::Proof,
        expression: &Expression<L::LoadedScalar>,
        evals: &HashMap<Query, L::LoadedScalar>,
        challenges: &[L::LoadedScalar],
        ys: &[&[L::LoadedScalar]],
        sum: &L::LoadedScalar,
        num_vars: usize,
        degree: usize,
    ) -> Result<Self::Output, Error> {
        let degree = expression.degree();
        let mut round_sum = sum.clone();
        let x = Vec::with_capacity(num_vars);
        for i in 0..num_vars {
            round_sum = SCR::verify(&proof.round_proofs[i], &round_sum, degree, i)?;
            x.push(proof.round_proofs[i].challenge);
        }

        let loader = sum.loader();
        loader.assert_eq(
            round_sum,
            evaluate::<C, L>(expression, num_vars, evals, challenges, ys, &x),
        );
        Ok(x)
    }
}

pub struct ClassicSumcheckProof<C, L> {
    round_proofs: Vec<ClassicSumcheckRoundProof<C, L>>,
}
