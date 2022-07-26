//! Implementing Structured Reference Strings for univariate polynomial KZG

use crate::{prelude::PCSErrors, StructuredReferenceString};
use ark_ec::{msm::FixedBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{end_timer, rand::RngCore, start_timer, One, UniformRand};

/// `UniversalParams` are the universal parameters for the KZG10 scheme.
/// Adapted from
/// https://github.com/arkworks-rs/poly-commit/blob/master/src/kzg10/data_structures.rs#L20
#[derive(Clone, Debug)]
pub struct UnivariateUniversalParams<E: PairingEngine> {
    /// Group elements of the form `{ \beta^i G }`, where `i` ranges from 0 to
    /// `degree`.
    pub powers_of_g: Vec<E::G1Affine>,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
}

/// Prover Parameters
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct UnivariateProverParam<C: AffineCurve> {
    pub powers_of_g: Vec<C>,
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
/// https://github.com/arkworks-rs/poly-commit/blob/master/src/kzg10/data_structures.rs#L236
#[derive(Clone, Debug)]
pub struct UnivariateVerifierParam<E: PairingEngine> {
    /// The generator of G1.
    pub g: E::G1Affine,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
}

impl<E: PairingEngine> StructuredReferenceString<E> for UnivariateUniversalParams<E> {
    type ProverParam = UnivariateProverParam<E::G1Affine>;
    type VerifierParam = UnivariateVerifierParam<E>;

    /// Extract the prover parameters from the public parameters.
    fn extract_prover_param(&self, supported_log_size: usize) -> Self::ProverParam {
        let support_size = 1usize << supported_log_size;
        let powers_of_g = self.powers_of_g[..=support_size].to_vec();

        Self::ProverParam { powers_of_g }
    }

    /// Extract the verifier parameters from the public parameters.
    fn extract_verifier_param(&self, _supported_log_size: usize) -> Self::VerifierParam {
        Self::VerifierParam {
            g: self.powers_of_g[0],
            h: self.h,
            beta_h: self.beta_h,
        }
    }

    /// Trim the universal parameters to specialize the public parameters
    /// for univariate polynomials to the given `supported_log_size`, and
    /// returns committer key and verifier key. `supported_log_size` should
    /// be in range `1..log(params.len())`
    fn trim(
        &self,
        supported_log_size: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSErrors> {
        let support_size = 1usize << supported_log_size;
        let powers_of_g = self.powers_of_g[..=support_size].to_vec();

        let pk = Self::ProverParam { powers_of_g };
        let vk = Self::VerifierParam {
            g: self.powers_of_g[0],
            h: self.h,
            beta_h: self.beta_h,
        };
        Ok((pk, vk))
    }

    /// Build SRS for testing.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: RngCore>(rng: &mut R, log_degree: usize) -> Result<Self, PCSErrors> {
        let max_degree = 1usize << log_degree;
        let setup_time = start_timer!(|| format!("KZG10::Setup with degree {}", max_degree));
        let beta = E::Fr::rand(rng);
        let g = E::G1Projective::rand(rng);
        let h = E::G2Projective::rand(rng);

        let mut powers_of_beta = vec![E::Fr::one()];

        let mut cur = beta;
        for _ in 0..max_degree {
            powers_of_beta.push(cur);
            cur *= &beta;
        }

        let window_size = FixedBaseMSM::get_mul_window_size(max_degree + 1);

        let scalar_bits = E::Fr::size_in_bits();
        let g_time = start_timer!(|| "Generating powers of G");
        // TODO: parallelization
        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g);
        let powers_of_g = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
            scalar_bits,
            window_size,
            &g_table,
            &powers_of_beta,
        );
        end_timer!(g_time);

        let powers_of_g = E::G1Projective::batch_normalization_into_affine(&powers_of_g);

        let h = h.into_affine();
        let beta_h = h.mul(beta).into_affine();

        let pp = Self {
            powers_of_g,
            h,
            beta_h,
        };
        end_timer!(setup_time);
        Ok(pp)
    }
}