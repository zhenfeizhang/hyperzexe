use crate::halo2_verifier::{loader::Loader, system::halo2::transcript};

#[allow(non_snake_case)]
#[derive(Clone, Debug)]
struct DotProductProof<L> {
    blind_eval: L::LoadedScalar,
    Qd: L::LoadedEcPoint,
    Ghat_d: L::LoadedEcPoint,
    z1: L::LoadedScalar,
    z2: L::LoadedScalar,
    bulletproof: BulletReductionProof<L>,
    c: L::LoadedScalar,
}

impl DotProductProof<L: Loader> {
    #[allow(non_snake_case)]
    fn read<T: TranscriptRead<C, L>>(n: usize, transcript: &mut T) -> Self {
        bulletproof = BulletReductionProof::read(n, transcript);
        let Qd = transcript.read_ec_point().unwrap();
        let Ghat_d = transcript.read_ec_point().unwrap();
        let c = transcript.squeeze_challenge();
        let zs = transcript.read_n_scalars(2).unwrap();

        Self {
            bulletproof,
            Qd,
            Ghat_d,
            c,
            z1: zs[0],
            z2: zs[1],
        }
    }

    fn verify(
        &self,
        gens: &DotProductProofGens<C>,
        b: &[L::LoadedScalar],
        Ga: &L::LoadedEcPoint,
        Qy: &L::LoadedEcPoint,
    ) -> Result<(), Error> {
        let n = gens.n;
        let loader = b[0].loader();
        assert_eq!(b.len(), n);

        let Gamma = *Ga + *Qy;

        let (g_hat, Gamma_hat, b_hat) = self.bullet_reduction_proof.verify(
            n,
            b,
            transcript,
            &Gamma.to_affine(),
            &gens.gens_n.G,
        )?;

        let lhs = (Gamma_hat.mul(c) + self.Ghat_d).to_affine().mul(b_hat) + self.Qd;
        let rhs = (g_hat.mul(self.z1) + gens.gens_1.G[0].mul(b_hat * self.z1))
            + gens.gens_1.h.mul(self.z2);

        loader.ec_point_assert_eq(lhs, rhs)
    }
}
