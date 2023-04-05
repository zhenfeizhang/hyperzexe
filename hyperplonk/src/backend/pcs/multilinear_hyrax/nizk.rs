#![allow(clippy::too_many_arguments)]

use super::{commitments::MultiCommitGens, math::Math, random::RandomTape};
use crate::{
    backend::{
        pcs::multilinear_hyrax::{commit_array, commit_element},
        util::transcript::TranscriptWrite,
    },
    Error,
};

mod bullet;
use bullet::BulletReductionProof;
use halo2_curves::{
    group::{ff::Field, Curve},
    CurveAffine,
};

#[derive(Debug, Clone)]
pub struct DotProductProofGens<C: CurveAffine> {
    pub n: usize,
    pub gens_n: MultiCommitGens<C>,
    pub gens_1: MultiCommitGens<C>,
}

impl<C: CurveAffine> DotProductProofGens<C> {
    pub fn new(n: usize, label: &[u8]) -> Result<Self, Error> {
        let (gens_n, gens_1) = MultiCommitGens::new(n + 1, label)?.split_at(n);
        Ok(DotProductProofGens { n, gens_n, gens_1 })
    }

    pub fn trim(&self, new_n: usize) -> Result<Self, Error> {
        if new_n > self.n {
            return Err(Error::InvalidPcsParam(format!(
                "SRS does not support target number of vars {}",
                new_n
            )));
        }

        let (gens_new_n, _) = self.gens_n.split_at(new_n);
        Ok(DotProductProofGens {
            n: new_n,
            gens_n: gens_new_n,
            gens_1: self.gens_1.clone(),
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DotProductProofLog<C: CurveAffine> {
    bullet_reduction_proof: BulletReductionProof<C>,
    Qd: C,
    Ghat_d: C,
    z1: C::Scalar,
    z2: C::Scalar,
}

impl<C: CurveAffine> DotProductProofLog<C> {
    fn protocol_name() -> &'static [u8] {
        b"dot product proof (log)"
    }

    pub fn prove(
        gens: &DotProductProofGens<C>,
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
        random_tape: &mut RandomTape<C::Scalar>,
        a_vec: &[C::Scalar],
        blind_a: &C::Scalar,
        b_vec: &[C::Scalar],
        y: &C::Scalar,
        blind_y: &C::Scalar,
        is_blind: bool,
    ) -> Result<(Self, C, C), Error> {
        let n = a_vec.len();
        assert_eq!(a_vec.len(), b_vec.len());
        assert_eq!(gens.n, n);
        assert_eq!(gens.gens_n.h, gens.gens_1.h);

        // produce randomness for generating a proof
        let (d, blind_Qd, blind_Ghat_d, blinds_vec) = if is_blind {
            let d = random_tape.random_scalar(b"d");
            let r_delta = random_tape.random_scalar(b"r_delta");
            let r_beta = random_tape.random_scalar(b"r_beta");

            let v1 = random_tape.random_vector(b"blinds_vec_1", 2 * n.log_2());
            let v2 = random_tape.random_vector(b"blinds_vec_2", 2 * n.log_2());
            let blinds_vec: Vec<(C::Scalar, C::Scalar)> =
                (0..v1.len()).map(|i| (v1[i], v2[i])).collect();

            (d, r_delta, r_beta, blinds_vec)
        } else {
            let v1 = vec![C::Scalar::zero(); 2 * n.log_2()];
            let v2 = vec![C::Scalar::zero(); 2 * n.log_2()];
            let blinds_vec = (0..v1.len())
                .map(|i| (v1[i], v2[i]))
                .collect::<Vec<(C::Scalar, C::Scalar)>>();
            (
                C::Scalar::zero(),
                C::Scalar::zero(),
                C::Scalar::zero(),
                blinds_vec,
            )
        };

        let blind_Gamma = *blind_a + *blind_y;

        let Ga = commit_array(a_vec, blind_a, &gens.gens_n).to_affine();
        let Qy = commit_element(y, blind_y, &gens.gens_1).to_affine();

        let (bullet_reduction_proof, _Gamma_hat, a_hat, b_hat, g_hat, blind_hat) =
            BulletReductionProof::<C>::prove(
                transcript,
                &gens.gens_1.G[0],
                &gens.gens_n.G,
                &gens.gens_n.h,
                a_vec,
                b_vec,
                &blind_Gamma,
                &blinds_vec,
            )?;
        let y_hat = a_hat * b_hat;

        let Qd = {
            let gens_hat = MultiCommitGens {
                n: 1,
                G: vec![g_hat],
                h: gens.gens_1.h,
            };
            commit_element(&d, &blind_Qd, &gens_hat).to_affine()
        };
        transcript.write_commitment(&Qd)?;

        let Ghat_d = commit_element(&d, &blind_Ghat_d, &gens.gens_1).to_affine();
        transcript.write_commitment(&Ghat_d)?;

        let c = transcript.squeeze_challenge();

        // Gamma_hat = a_hat[G_hat] + y_hat[Q] + blind_hat[H]
        // b_hat * Gamma_hat = y_hat[G_hat] + b * y_hat[Q] + b * blind_hat[H]
        // The goal is to prove that there exists y_hat that satisfies the above
        // equation. This can be proved by ZKP for discrete-log identity.
        let z1 = d + c * y_hat;
        let z2 = b_hat * (c * blind_hat + blind_Ghat_d) + blind_Qd;
        Ok((
            DotProductProofLog {
                bullet_reduction_proof,
                Qd,
                Ghat_d,
                z1,
                z2,
            },
            Ga,
            Qy,
        ))
    }

    pub fn verify(
        &self,
        n: usize,
        gens: &DotProductProofGens<C>,
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
        b: &[C::Scalar],
        Ga: &C,
        Qy: &C,
    ) -> Result<bool, Error> {
        assert_eq!(gens.n, n);
        assert_eq!(b.len(), n);

        let Gamma = *Ga + *Qy;

        let (g_hat, Gamma_hat, b_hat) = self.bullet_reduction_proof.verify(
            n,
            b,
            transcript,
            &Gamma.to_affine(),
            &gens.gens_n.G,
        )?;
        transcript.write_commitment(&self.Qd)?;
        transcript.write_commitment(&self.Ghat_d)?;

        let c = transcript.squeeze_challenge();

        let lhs = (Gamma_hat.mul(c) + self.Ghat_d).to_affine().mul(b_hat) + self.Qd;
        let rhs = (g_hat.mul(self.z1) + gens.gens_1.G[0].mul(b_hat * self.z1))
            + gens.gens_1.h.mul(self.z2);

        Ok(lhs == rhs)
    }
}

#[cfg(test)]
mod tests {
    use crate::backend::util::{arithmetic::sum, test::std_rng, transcript::Keccak256Transcript};

    use super::*;
    use halo2_curves::bn256::{Fr, G1Affine as C};

    #[test]
    fn check_dotproductproof_log() -> Result<(), Error> {
        let mut csprng = std_rng();

        let n = 1024;

        let gens = DotProductProofGens::<C>::new(n, b"test-1024")?;

        let a: Vec<Fr> = (0..n).map(|_i| Fr::random(&mut csprng)).collect();
        let b: Vec<Fr> = (0..n).map(|_i| Fr::random(&mut csprng)).collect();
        let y: Fr = sum((0..a.len()).map(|i| a[i] * b[i]));

        let r_a = Fr::random(&mut csprng);
        let r_y = Fr::random(&mut csprng);

        let mut random_tape = RandomTape::new(b"proof");
        let mut prover_transcript = Keccak256Transcript::<Vec<u8>>::default();

        let (proof, Ca, Cy) = DotProductProofLog::prove(
            &gens,
            &mut prover_transcript,
            &mut random_tape,
            &a,
            &r_a,
            &b,
            &y,
            &r_y,
            true,
        )?;
        let mut verifier_transcript = Keccak256Transcript::<Vec<u8>>::default();
        assert!(proof.verify(n, &gens, &mut verifier_transcript, &b, &Ca, &Cy)?);
        Ok(())
    }
}
