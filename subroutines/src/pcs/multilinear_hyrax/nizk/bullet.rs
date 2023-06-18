//! This module is an adaptation of code from the bulletproofs crate.
//! See NOTICE.md for more details
#![allow(non_snake_case)]
#![allow(clippy::type_complexity)]
#![allow(clippy::too_many_arguments)]
use crate::PCSError;

use super::super::math::Math;
use ark_ec::{msm::VariableBaseMSM, AffineCurve, ProjectiveCurve};
use ark_ff::{batch_inversion, Field, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{One, Zero};
use core::iter;
use transcript::IOPTranscript;

#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize, PartialEq, Eq)]
pub struct BulletReductionProof<C: AffineCurve> {
    pub L_vec: Vec<C>,
    pub R_vec: Vec<C>,
}

impl<C: AffineCurve> BulletReductionProof<C> {
    /// Create an inner-product proof.
    ///
    /// The proof is created with respect to the bases \\(G\\).
    ///
    /// The `transcript` is passed in as a parameter so that the
    /// challenges depend on the *entire* transcript (including parent
    /// protocols).
    ///
    /// The lengths of the vectors must all be the same, and must all be
    /// either 0 or a power of 2.
    pub fn prove(
        transcript: &mut IOPTranscript<C::ScalarField>,
        Q: &C,
        G_vec: &[C],
        H: &C,
        a_vec: &[C::ScalarField],
        b_vec: &[C::ScalarField],
        blind: &C::ScalarField,
        blinds_vec: &[(C::ScalarField, C::ScalarField)],
    ) -> Result<(Self, C, C::ScalarField, C::ScalarField, C, C::ScalarField), PCSError> {
        // Create slices G, H, a, b backed by their respective
        // vectors.  This lets us reslice as we compress the lengths
        // of the vectors in the main loop below.
        let mut G: &mut [C] = &mut G_vec.to_owned()[..];
        let mut a: &mut [C::ScalarField] = &mut a_vec.to_owned()[..];
        let mut b: &mut [C::ScalarField] = &mut b_vec.to_owned()[..];

        // All of the input vectors must have a length that is a power of two.
        let mut n = G.len();
        assert!(n.is_power_of_two());
        let lg_n = n.log_2();

        // All of the input vectors must have the same length.
        assert_eq!(G.len(), n);
        assert_eq!(a.len(), n);
        assert_eq!(b.len(), n);
        assert_eq!(blinds_vec.len(), 2 * lg_n);

        let mut L_vec = Vec::with_capacity(lg_n);
        let mut R_vec = Vec::with_capacity(lg_n);
        let mut blinds_iter = blinds_vec.iter();
        let mut blind_fin = *blind;

        while n != 1 {
            n /= 2;
            let (a_L, a_R): (&mut [C::ScalarField], &mut [C::ScalarField]) = a.split_at_mut(n);
            let (b_L, b_R): (&mut [C::ScalarField], &mut [C::ScalarField]) = b.split_at_mut(n);
            let (G_L, G_R) = G.split_at_mut(n);

            let c_L = inner_product(a_L, b_R);
            let c_R = inner_product(a_R, b_L);

            let (blind_L, blind_R) = blinds_iter.next().unwrap();

            let L = {
                let a_L_bigint: Vec<_> = a_L.iter().map(|x| x.into_repr()).collect();
                (VariableBaseMSM::multi_scalar_mul(G_R, a_L_bigint.as_slice())
                    + Q.mul(c_L)
                    + H.mul(*blind_L))
                .into_affine()
            };

            let R = {
                let a_R_bigint: Vec<_> = a_R.iter().map(|x| x.into_repr()).collect();
                (VariableBaseMSM::multi_scalar_mul(G_L, a_R_bigint.as_slice())
                    + Q.mul(c_R)
                    + H.mul(*blind_R))
                .into_affine()
            };

            transcript.append_serializable_element(b"L", &L)?;
            transcript.append_serializable_element(b"R", &R)?;

            let u = transcript.get_and_append_challenge(b"u")?;
            assert!(u != C::ScalarField::zero());
            let u_inv = u.inverse().unwrap();
            {
                let G_L_vec = {
                    let mut G_L_projective = vec![C::zero().into_projective(); n];
                    for i in 0..n {
                        a_L[i] = a_L[i] * u + u_inv * a_R[i];
                        b_L[i] = b_L[i] * u_inv + u * b_R[i];
                        G_L_projective[i] = G_L[i].mul(u_inv) + G_R[i].mul(u);
                    }
                    C::Projective::batch_normalization_into_affine(&G_L_projective)
                };
                for i in 0..n {
                    G_L[i] = G_L_vec[i];
                    assert!(G_L[i].is_zero() == false);
                }
            }

            blind_fin = blind_fin + *blind_L * u * u + *blind_R * u_inv * u_inv;

            L_vec.push(L);
            R_vec.push(R);

            a = a_L;
            b = b_L;
            G = G_L;
        }

        let Gamma_hat = VariableBaseMSM::multi_scalar_mul(
            &[G[0], *Q, *H],
            &[
                a[0].into_repr(),
                (a[0] * b[0]).into_repr(),
                blind_fin.into_repr(),
            ],
        )
        .into_affine();

        Ok((
            BulletReductionProof { L_vec, R_vec },
            Gamma_hat,
            a[0],
            b[0],
            G[0],
            blind_fin,
        ))
    }

    /// Computes three vectors of verification scalars \\([u\_{i}^{2}]\\),
    /// \\([u\_{i}^{-2}]\\) and \\([s\_{i}]\\) for combined multiscalar
    /// multiplication in a parent protocol. See [inner product protocol
    /// notes](index.html#verification-equation) for details. The verifier
    /// must provide the input length \\(n\\) explicitly to avoid unbounded
    /// allocation within the inner product proof.
    fn verification_scalars(
        &self,
        n: usize,
        transcript: &mut IOPTranscript<C::ScalarField>,
    ) -> Result<
        (
            Vec<C::ScalarField>,
            Vec<C::ScalarField>,
            Vec<C::ScalarField>,
        ),
        PCSError,
    > {
        let lg_n = self.L_vec.len();
        if lg_n >= 32 {
            // 4 billion multiplications should be enough for anyone
            // and this check prevents overflow in 1<<lg_n below.
            return Err(PCSError::InvalidVerifier(
                "Inner product proof is too large".to_string(),
            ));
        }
        if n != (1 << lg_n) {
            return Err(PCSError::InvalidVerifier(format!(
                "Inner product proof is for {} elements, but verifier provided {}",
                1 << lg_n,
                n
            )));
        }

        // 1. Recompute x_k,...,x_1 based on the proof transcript
        let mut challenges = Vec::with_capacity(lg_n);
        for (L, R) in self.L_vec.iter().zip(self.R_vec.iter()) {
            transcript.append_serializable_element(b"L", L)?;
            transcript.append_serializable_element(b"R", R)?;
            challenges.push(transcript.get_and_append_challenge(b"u")?);
        }

        // 2. Compute 1/(u_k...u_1) and 1/u_k, ..., 1/u_1
        let mut challenges_inv = challenges.clone();
        batch_inversion(&mut challenges_inv);
        let allinv = challenges_inv
            .iter()
            .fold(C::ScalarField::one(), |acc, x| acc * x);

        // 3. Compute u_i^2 and (1/u_i)^2
        for i in 0..lg_n {
            challenges[i] = challenges[i].square();
            challenges_inv[i] = challenges_inv[i].square();
        }
        let challenges_sq = challenges;
        let challenges_inv_sq = challenges_inv;

        // 4. Compute s values inductively.
        let mut s = Vec::with_capacity(n);
        s.push(allinv);
        for i in 1..n {
            let lg_i = (32 - 1 - (i as u32).leading_zeros()) as usize;
            let k = 1 << lg_i;
            // The challenges are stored in "creation order" as [u_k,...,u_1],
            // so u_{lg(i)+1} = is indexed by (lg_n-1) - lg_i
            let u_lg_i_sq = challenges_sq[(lg_n - 1) - lg_i];
            s.push(s[i - k] * u_lg_i_sq);
        }

        Ok((challenges_sq, challenges_inv_sq, s))
    }

    /// This method is for testing that proof generation work,
    /// but for efficiency the actual protocols would use `verification_scalars`
    /// method to combine inner product verification with other checks
    /// in a single multiscalar multiplication.
    pub fn verify(
        &self,
        n: usize,
        a: &[C::ScalarField],
        transcript: &mut IOPTranscript<C::ScalarField>,
        Gamma: &C,
        G: &[C],
    ) -> Result<(C, C, C::ScalarField), PCSError> {
        let (u_sq, u_inv_sq, s) = self.verification_scalars(n, transcript)?;

        let a_hat = inner_product(a, s.as_slice());
        let G_hat = VariableBaseMSM::multi_scalar_mul(
            G,
            s.into_iter()
                .map(|x| x.into_repr())
                .collect::<Vec<_>>()
                .as_slice(),
        )
        .into_affine();

        let points: Vec<C> = self
            .L_vec
            .iter()
            .chain(self.R_vec.iter())
            .chain(iter::once(Gamma))
            .cloned()
            .collect();

        let scalars: Vec<_> = u_sq
            .iter()
            .chain(u_inv_sq.iter())
            .chain(iter::once(&C::ScalarField::one()))
            .map(|x| x.into_repr())
            .collect();

        let Gamma_hat =
            VariableBaseMSM::multi_scalar_mul(points.as_slice(), scalars.as_slice()).into_affine();

        Ok((G_hat, Gamma_hat, a_hat))
    }
}

/// Computes an inner product of two vectors
/// \\[
///    {\langle {\mathbf{a}}, {\mathbf{b}} \rangle} = \sum\_{i=0}^{n-1} a\_i
/// \cdot b\_i. \\]
/// Panics if the lengths of \\(\mathbf{a}\\) and \\(\mathbf{b}\\) are not
/// equal.
pub fn inner_product<F: Field>(a: &[F], b: &[F]) -> F {
    assert!(
        a.len() == b.len(),
        "inner_product(a,b): lengths of vectors do not match"
    );
    let mut out = F::zero();
    for i in 0..a.len() {
        out += a[i] * b[i];
    }
    out
}
