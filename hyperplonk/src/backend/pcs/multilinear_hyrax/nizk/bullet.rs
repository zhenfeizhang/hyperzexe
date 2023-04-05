//! This module is an adaptation of code from the bulletproofs crate.
//! See NOTICE.md for more details
#![allow(non_snake_case)]
#![allow(clippy::type_complexity)]
#![allow(clippy::too_many_arguments)]

use halo2_curves::{
    group::{
        ff::{BatchInvert, Field},
        Curve,
    },
    CurveAffine,
};

use crate::{
    backend::util::{arithmetic::variable_base_msm, transcript::TranscriptWrite},
    Error,
};

use super::super::math::Math;
use core::iter;
use std::iter::once;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BulletReductionProof<C: CurveAffine> {
    pub L_vec: Vec<C>,
    pub R_vec: Vec<C>,
}

impl<C: CurveAffine> BulletReductionProof<C> {
    /// Create an inner-product proof.
    ///
    /// If <a, b> = y, this is to prove the folding scheme runs correctly to
    /// prove Gamma = <a, [G]> + <y, [Q]>. Where [G] is a base vector to
    /// commit n-sized vectors and [Q] is a base element to commit a single
    /// element. [H] is the base for the blinding factor.
    ///
    /// At the final step
    /// Gamma_hat = sum_{1 <= i <= n}( L_i * u_i^2 + R_i * u_i^{-2} ) + <a, [G]>
    /// + blind_a[H] + <y, [Q]> + blind_y[H]           = sum_{1 <= i <= n}(
    /// (<a_Li, [G_Ri]> + <a_Li, b_Ri>[Q] + blind_Li[H]) * u_i^2 )
    ///           + sum_{1 <= i <= n}( (<a_Ri, [G_Li]> + <a_Ri, b_Li>[Q] +
    /// blind_Ri[H]) * u_i^{-2} )           + <a, [G]> + blind_a[H] + <y,
    /// [Q]> + blind_y[H]           = <a_hat, [G_hat]> + (a_hat * b_hat)[Q]
    /// + blind_hat[H]
    ///
    /// where a_hat = sum_{0 <= i < 2^n} a_i * prod_{0 <= j < n} u_i^{
    /// (-1)^{i[j]} }       b_hat = sum_{0 <= i < 2^n} b_i * prod_{0 <= j <
    /// n} u_i^{ (-1)^{1 - i[j]} }       G_hat = sum_{0 <= i < 2^n} G_i *
    /// prod_{0 <= j < n} u_i^{ (-1)^{1 - i[j]} }       blind_hat =
    /// blind_a[H] + blind_y[H] + sum_{1 <= i <= n}( blind_Li[H] * u^2 +
    /// blind_Ri[H] * u^{-2} )
    ///
    /// This function will return Gamma_hat, a_hat, b_hat, G_hat, blind_hat
    pub fn prove(
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
        Q: &C,
        G_vec: &[C],
        H: &C,
        a_vec: &[C::Scalar],
        b_vec: &[C::Scalar],
        blind: &C::Scalar,                     // blind = blind_a + blind_y
        blinds_vec: &[(C::Scalar, C::Scalar)], // (blind_Li, blind_Ri)_{1 <= i <= n}
    ) -> Result<(Self, C, C::Scalar, C::Scalar, C, C::Scalar), Error> {
        // Create slices G, H, a, b backed by their respective
        // vectors.  This lets us reslice as we compress the lengths
        // of the vectors in the main loop below.
        let mut G: &mut [C] = &mut G_vec.to_owned()[..];
        let mut a: &mut [C::Scalar] = &mut a_vec.to_owned()[..];
        let mut b: &mut [C::Scalar] = &mut b_vec.to_owned()[..];

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
        let mut blind_hat = *blind;

        while n != 1 {
            n /= 2;
            let (a_L, a_R): (&mut [C::Scalar], &mut [C::Scalar]) = a.split_at_mut(n);
            let (b_L, b_R): (&mut [C::Scalar], &mut [C::Scalar]) = b.split_at_mut(n);
            let (G_L, G_R) = G.split_at_mut(n);

            let c_L = inner_product(a_L, b_R);
            let c_R = inner_product(a_R, b_L);

            let (blind_L, blind_R) = blinds_iter.next().unwrap();

            // L = <a_L, [G_R]> + <a_L, b_R>[Q] + blind_L[H]
            let L = {
                let scalars = a_L.iter().chain(once(&c_L)).chain(once(blind_L));
                let points = G_R.iter().chain(once(Q)).chain(once(H));
                variable_base_msm(scalars, points).to_affine()
            };

            // R = <a_R, [G_L]> + <a_R, b_L>[Q] + blind_R[H]
            let R = {
                let scalars = a_R.iter().chain(once(&c_R)).chain(once(blind_R));
                let points = G_L.iter().chain(once(Q)).chain(once(H));
                variable_base_msm(scalars, points).to_affine()
            };

            transcript.write_commitment(&L)?;
            transcript.write_commitment(&R)?;

            // Update   a = a_L * u + a_R * u_inv,
            //          b = b_L * u_inv + a_R * u,
            //          G = G_L * u_inv + G_R * u
            let u = transcript.squeeze_challenge();
            assert!(u != C::Scalar::zero());
            let u_inv = u.invert().unwrap();
            {
                let mut G_L_vec: Vec<C> = vec![C::identity(); n];
                let mut G_L_projective: Vec<C::CurveExt> = Vec::with_capacity(n);
                for i in 0..n {
                    a_L[i] = a_L[i] * u + u_inv * a_R[i];
                    b_L[i] = b_L[i] * u_inv + u * b_R[i];
                    G_L_projective.push(G_L[i].mul(u_inv) + G_R[i].mul(u));
                }
                C::Curve::batch_normalize(&G_L_projective, &mut G_L_vec);
                for i in 0..n {
                    G_L[i] = G_L_vec[i];
                    assert!(!bool::from(G_L[i].is_identity()));
                }
            }

            blind_hat = blind_hat + *blind_L * u * u + *blind_R * u_inv * u_inv;

            L_vec.push(L);
            R_vec.push(R);

            a = a_L;
            b = b_L;
            G = G_L;
        }

        let Gamma_hat =
            variable_base_msm(&[a[0], (a[0] * b[0]), blind_hat], &[G[0], *Q, *H]).to_affine();

        Ok((
            BulletReductionProof { L_vec, R_vec },
            Gamma_hat,
            a[0],
            b[0],
            G[0],
            blind_hat,
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
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
    ) -> Result<(Vec<C::Scalar>, Vec<C::Scalar>, Vec<C::Scalar>), Error> {
        let lg_n = self.L_vec.len();
        if lg_n >= 32 {
            // 4 billion multiplications should be enough for anyone
            // and this check prevents overflow in 1<<lg_n below.
            return Err(Error::InvalidPcsOpen(
                "Inner product proof is too large".to_string(),
            ));
        }
        if n != (1 << lg_n) {
            return Err(Error::InvalidPcsOpen(format!(
                "Inner product proof is for {} elements, but verifier provided {}",
                1 << lg_n,
                n
            )));
        }

        // 1. Recompute x_k,...,x_1 based on the proof transcript
        let mut challenges = Vec::with_capacity(lg_n);
        for (L, R) in self.L_vec.iter().zip(self.R_vec.iter()) {
            transcript.write_commitment(L)?;
            transcript.write_commitment(R)?;
            challenges.push(transcript.squeeze_challenge());
        }

        // 2. Compute 1/(u_k...u_1) and 1/u_k, ..., 1/u_1
        let mut challenges_inv = challenges.clone();
        challenges_inv.iter_mut().batch_invert();
        let allinv = challenges_inv
            .iter()
            .fold(C::Scalar::one(), |acc, x| acc * x);

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
        a: &[C::Scalar],
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
        Gamma: &C,
        G: &[C],
    ) -> Result<(C, C, C::Scalar), Error> {
        let (u_sq, u_inv_sq, s) = self.verification_scalars(n, transcript)?;

        let a_hat = inner_product(a, s.as_slice());
        let G_hat = variable_base_msm(&s, G).to_affine();

        let points = self
            .L_vec
            .iter()
            .chain(self.R_vec.iter())
            .chain(once(Gamma));

        let one = C::Scalar::one();
        let scalars = u_sq.iter().chain(u_inv_sq.iter()).chain(iter::once(&one));

        let Gamma_hat = variable_base_msm(scalars, points).to_affine();

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
