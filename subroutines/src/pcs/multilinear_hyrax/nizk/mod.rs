#![allow(clippy::too_many_arguments)]

use super::{commitments::MultiCommitGens, math::Math, random::RandomTape};
use crate::pcs::multilinear_hyrax::commit_array;
use crate::{pcs::multilinear_hyrax::commit_element, PCSError};
use ark_ec::AffineCurve;
use ark_ec::ProjectiveCurve;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::Zero;
use transcript::IOPTranscript;

mod bullet;
use bullet::BulletReductionProof;

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct DotProductProofGens<C: AffineCurve> {
    n: usize,
    pub gens_n: MultiCommitGens<C>,
    pub gens_1: MultiCommitGens<C>,
}

impl<C: AffineCurve> DotProductProofGens<C> {
    pub fn new(n: usize, label: &[u8]) -> Result<Self, PCSError> {
        let (gens_n, gens_1) = MultiCommitGens::new(n + 1, label)?.split_at(n);
        Ok(DotProductProofGens { n, gens_n, gens_1 })
    }

    pub fn trim(&self, new_n: usize) -> Result<Self, PCSError> {
        if new_n > self.n {
            return Err(PCSError::InvalidParameters(format!(
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

#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize, PartialEq, Eq)]
pub struct DotProductProofLog<C: AffineCurve> {
    bullet_reduction_proof: BulletReductionProof<C>,
    delta: C,
    beta: C,
    z1: C::ScalarField,
    z2: C::ScalarField,
}

impl<C: AffineCurve> DotProductProofLog<C> {
    fn protocol_name() -> &'static [u8] {
        b"dot product proof (log)"
    }

    pub fn prove(
        gens: &DotProductProofGens<C>,
        transcript: &mut IOPTranscript<C::ScalarField>,
        random_tape: &mut RandomTape<C::ScalarField>,
        x_vec: &[C::ScalarField],
        blind_x: &C::ScalarField,
        a_vec: &[C::ScalarField],
        y: &C::ScalarField,
        blind_y: &C::ScalarField,
        is_blind: bool,
    ) -> Result<(Self, C, C), PCSError> {
        transcript.append_protocol_name(DotProductProofLog::<C>::protocol_name())?;

        let n = x_vec.len();
        assert_eq!(x_vec.len(), a_vec.len());
        assert_eq!(gens.n, n);
        assert_eq!(gens.gens_n.h, gens.gens_1.h);

        // produce randomness for generating a proof
        // TODO: double check whether this is correct.
        let (d, r_delta, r_beta, blinds_vec) = if is_blind {
            let d = random_tape.random_scalar(b"d");
            let r_delta = random_tape.random_scalar(b"r_delta");
            let r_beta = random_tape.random_scalar(b"r_beta");

            let v1 = random_tape.random_vector(b"blinds_vec_1", 2 * n.log_2());
            let v2 = random_tape.random_vector(b"blinds_vec_2", 2 * n.log_2());
            let blinds_vec = (0..v1.len())
                .map(|i| (v1[i], v2[i]))
                .collect::<Vec<(C::ScalarField, C::ScalarField)>>();

            (d, r_delta, r_beta, blinds_vec)
        } else {
            let v1 = vec![C::ScalarField::zero(); 2 * n.log_2()];
            let v2 = vec![C::ScalarField::zero(); 2 * n.log_2()];
            let blinds_vec = (0..v1.len())
                .map(|i| (v1[i], v2[i]))
                .collect::<Vec<(C::ScalarField, C::ScalarField)>>();
            (
                C::ScalarField::zero(),
                C::ScalarField::zero(),
                C::ScalarField::zero(),
                blinds_vec,
            )
        };

        let blind_Gamma = *blind_x + *blind_y;

        let Cx = commit_array(x_vec, blind_x, &gens.gens_n).into_affine();
        transcript.append_serializable_element(b"Cx", &Cx)?;
        let Cy = commit_element(y, blind_y, &gens.gens_1).into_affine();
        transcript.append_serializable_element(b"Cy", &Cy)?;

        transcript.append_field_vectors(b"a", a_vec)?;

        let (bullet_reduction_proof, _Gamma_hat, x_hat, a_hat, g_hat, rhat_Gamma) =
            BulletReductionProof::<C>::prove(
                transcript,
                &gens.gens_1.G[0],
                &gens.gens_n.G,
                &gens.gens_n.h,
                x_vec,
                a_vec,
                &blind_Gamma,
                &blinds_vec,
            )?;
        let y_hat = x_hat * a_hat;

        let delta = {
            let gens_hat = MultiCommitGens {
                n: 1,
                G: vec![g_hat],
                h: gens.gens_1.h,
            };
            commit_element(&d, &r_delta, &gens_hat).into_affine()
        };
        transcript.append_serializable_element(b"delta", &delta)?;

        let beta = commit_element(&d, &r_beta, &gens.gens_1).into_affine();
        transcript.append_serializable_element(b"beta", &beta)?;

        let c = transcript.get_and_append_challenge(b"c")?;

        let z1 = d + c * y_hat;
        let z2 = a_hat * (c * rhat_Gamma + r_beta) + r_delta;
        Ok((
            DotProductProofLog {
                bullet_reduction_proof,
                delta,
                beta,
                z1,
                z2,
            },
            Cx,
            Cy,
        ))
    }

    pub fn verify(
        &self,
        n: usize,
        gens: &DotProductProofGens<C>,
        transcript: &mut IOPTranscript<C::ScalarField>,
        a: &[C::ScalarField],
        Cx: &C,
        Cy: &C,
    ) -> Result<bool, PCSError> {
        assert_eq!(gens.n, n);
        assert_eq!(a.len(), n);

        transcript.append_protocol_name(DotProductProofLog::<C>::protocol_name())?;
        transcript.append_serializable_element(b"Cx", Cx)?;
        transcript.append_serializable_element(b"Cy", Cy)?;
        transcript.append_field_vectors(b"a", a)?;

        let Gamma = *Cx + *Cy;

        let (g_hat, Gamma_hat, a_hat) =
            self.bullet_reduction_proof
                .verify(n, a, transcript, &Gamma, &gens.gens_n.G)?;
        transcript.append_serializable_element(b"delta", &self.delta)?;
        transcript.append_serializable_element(b"beta", &self.beta)?;

        let c = transcript.get_and_append_challenge(b"c")?;

        let lhs = (Gamma_hat.mul(c).into_affine() + self.beta)
            .mul(a_hat)
            .into_affine()
            + self.delta;
        let rhs = ((g_hat.mul(self.z1) + gens.gens_1.G[0].mul(a_hat * self.z1))
            + gens.gens_1.h.mul(self.z2))
        .into_affine();

        Ok(lhs == rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;
    use ark_bls12_381::G1Affine;
    use ark_ff::UniformRand;
    use ark_std::test_rng;

    type Scalar = Fr;
    type Curve = G1Affine;

    #[test]
    fn check_dotproductproof_log() -> Result<(), PCSError> {
        let mut csprng = test_rng();

        let n = 1024;

        let gens = DotProductProofGens::<Curve>::new(n, b"test-1024")?;

        let x: Vec<Scalar> = (0..n).map(|_i| Scalar::rand(&mut csprng)).collect();
        let a: Vec<Scalar> = (0..n).map(|_i| Scalar::rand(&mut csprng)).collect();
        let y: Scalar = (0..x.len()).map(|i| x[i] * a[i]).sum();

        let r_x = Scalar::rand(&mut csprng);
        let r_y = Scalar::rand(&mut csprng);

        let mut random_tape = RandomTape::new(b"proof");
        let mut prover_transcript = IOPTranscript::<Scalar>::new(b"example");

        let (proof, Cx, Cy) = DotProductProofLog::prove(
            &gens,
            &mut prover_transcript,
            &mut random_tape,
            &x,
            &r_x,
            &a,
            &y,
            &r_y,
            true,
        )?;
        let mut verifier_transcript = IOPTranscript::<Scalar>::new(b"example");
        assert!(proof.verify(n, &gens, &mut verifier_transcript, &a, &Cx, &Cy)?);
        Ok(())
    }
}
