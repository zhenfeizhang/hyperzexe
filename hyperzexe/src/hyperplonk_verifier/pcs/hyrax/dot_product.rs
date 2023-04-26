use halo2_proofs::curves::CurveAffine;
use hyperplonk::backend::pcs::prelude::DotProductProofGens;

use crate::{
    halo2_verifier::{loader::Loader, util::transcript::TranscriptRead},
    hyperplonk_verifier::Error,
};

use super::bulletproof::BulletReductionProof;

#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub(super) struct DotProductProof<C: CurveAffine, L: Loader<C>> {
    Qd: L::LoadedEcPoint,
    Ghat_d: L::LoadedEcPoint,
    z1: L::LoadedScalar,
    z2: L::LoadedScalar,
    bulletproof: BulletReductionProof<C, L>,
    c: L::LoadedScalar,
}

#[allow(non_snake_case)]
impl<C: CurveAffine, L: Loader<C>> DotProductProof<C, L> {
    pub(super) fn read_proof<T: TranscriptRead<C, L>>(n: usize, transcript: &mut T) -> Self {
        let bulletproof = BulletReductionProof::read(n, transcript);
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

    pub(super) fn verify(
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

        let (g_hat, Gamma_hat, b_hat) =
            self.bulletproof
                .verify(n, &b, &Gamma.to_affine(), &gens.gens_n.G)?;

        let lhs = (Gamma_hat * self.c + self.Ghat_d).to_affine().mul(b_hat) + self.Qd;
        let rhs = (g_hat.mul(self.z1) + gens.gens_1.G[0].mul(b_hat * self.z1))
            + gens.gens_1.h.mul(self.z2);

        loader.ec_point_assert_eq(lhs, rhs)
    }
}
