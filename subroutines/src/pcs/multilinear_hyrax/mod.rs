//! Main module for multilinear Hyrax commitment scheme (based on implementation
//! in Spartan)
#![allow(non_snake_case)]

use self::{
    batching::{batch_verify_single_point_internal, BatchProof},
    dense_mlpoly::{MultilinearHyraxProof, PolyCommitmentBlinds, PolyCommitmentGens}, nizk::DotProductProofGens,
};
use crate::{
    pcs::{
        multilinear_hyrax::{
            batching::multi_open_single_point_internal,
            commitments::{commit_array, commit_element, MultiCommitGens},
            math::Math,
            random::RandomTape,
        },
        structs::HyraxCommitment,
        PCSError,
    },
    PolynomialCommitmentScheme,
};

use arithmetic::evaluate_opt;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_poly::DenseMultilinearExtension;
use ark_std::{
    borrow::Borrow,
    end_timer,
    marker::PhantomData,
    rand::{CryptoRng, RngCore},
    start_timer,
    sync::Arc,
    vec,
    vec::Vec,
    Zero,
};
use transcript::IOPTranscript;

pub(crate) mod batching;
mod commitments;
mod dense_mlpoly;
mod math;
mod nizk;
mod random;

/// Hyrax Polynomial Commitment Scheme on multilinear polynomials.
pub struct MultilinearHyraxPCS<C: AffineCurve> {
    #[doc(hidden)]
    _curve: PhantomData<C>,
}

pub type HyraxSRS<C> = DotProductProofGens<C>;
pub type HyraxProverParam<C> = PolyCommitmentGens<C>;
pub type HyraxVerifierParam<C> = PolyCommitmentGens<C>;

impl<C> PolynomialCommitmentScheme<C> for MultilinearHyraxPCS<C>
where
    C: AffineCurve,
{
    // Parameters
    type ProverParam = HyraxProverParam<C>;
    type VerifierParam = HyraxVerifierParam<C>;
    /// Structured reference string
    type SRS = HyraxSRS<C>;
    // Polynomial and its associated types
    type Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>;
    type Point = Vec<C::ScalarField>;
    // Commitments and proofs
    type Commitment = HyraxCommitment<C>;
    type Proof = MultilinearHyraxProof<C>;
    type BatchProof = BatchProof<C, Self>;

    /// Build SRS for testing.
    ///
    /// - For univariate polynomials, `log_size` is the log of maximum degree.
    /// - For multilinear polynomials, `log_size` is the number of variables.
    ///
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: RngCore + CryptoRng>(
        _rng: &mut R,
        log_size: usize,
    ) -> Result<Self::SRS, PCSError> {
        let timer = start_timer!(|| "MultilinearHyraxPCS::gen_srs_for_testing");
        if log_size == 0 {
            return Err(PCSError::InvalidParameters(
                "constant polynomial not supported".to_string(),
            ));
        }
        let right_num_vars = log_size + 1 >> 1;
        let res = DotProductProofGens::new(1 << right_num_vars, b"new srs")?;
        end_timer!(timer);
        Ok(res)
    }

    /// Trim the universal parameters to specialize the public parameters.
    /// Input both `supported_log_degree` for univariate and
    /// `supported_num_vars` for multilinear.
    fn trim(
        srs: impl Borrow<Self::SRS>,
        supported_degree: Option<usize>,
        supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        let timer = start_timer!(|| "MultilinearHyraxPCS::trim");
        assert!(supported_degree.is_none());
        let supported_num_vars = match supported_num_vars {
            Some(p) => p,
            None => {
                return Err(PCSError::InvalidParameters(
                    "multilinear should receive a num_var param".to_string(),
                ))
            },
        };
        let param = PolyCommitmentGens::new_from_srs(srs.borrow(), supported_num_vars)?;
        end_timer!(timer);

        Ok((param.clone(), param))
    }

    /// Generate a commitment for a polynomial.
    ///
    /// This function takes `2^num_vars` number of scalar multiplications over
    /// G1.
    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, PCSError> {
        let prover_param = prover_param.borrow();
        let timer = start_timer!(|| "MultilinearHyraxPCS::commit");
        let n = poly.num_vars.pow2();
        let ell = poly.num_vars;

        let (left_num_vars, right_num_vars) = prover_param.compute_factored_lens(ell);
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();
        assert_eq!(L_size * R_size, n);

        // We don't need hiding property here.
        let blinds: PolyCommitmentBlinds<C> = PolyCommitmentBlinds {
            blinds: vec![C::ScalarField::zero(); L_size],
        };

        let res = Self::commit_inner(poly, &blinds.blinds, &prover_param.gens.gens_n);
        end_timer!(timer);
        Ok(res)
    }

    /// Generate a commitment for polynomials
    fn multi_commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &[Self::Polynomial],
    ) -> Result<Self::Commitment, PCSError> {
        let prover_param = prover_param.borrow();
        let timer = start_timer!(|| "MultilinearHyraxPCS::multi_commit");

        // Concatenate all polynomial evaluations and pad to the next power of 2.
        let mut poly_evals = Vec::new();
        for p in poly {
            poly_evals.extend(&p.evaluations);
        }
        poly_evals.extend(vec![
            C::ScalarField::zero();
            poly_evals.len().next_power_of_two() - poly_evals.len()
        ]);

        let num_vars = poly_evals.len().log_2();
        let poly = DenseMultilinearExtension::from_evaluations_vec(num_vars, poly_evals);

        let (left_num_vars, _) = prover_param.compute_factored_lens(num_vars);
        let L_size = left_num_vars.pow2();

        // We don't need hiding property here.
        let blinds: PolyCommitmentBlinds<C> = PolyCommitmentBlinds {
            blinds: vec![C::ScalarField::zero(); L_size],
        };

        let res = Self::commit_inner(&poly, &blinds.blinds, &prover_param.gens.gens_n);
        end_timer!(timer);
        Ok(res)
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same. This function does not need to take the evaluation value as an
    /// input.
    ///
    /// This function takes 2^{num_var +1} number of scalar multiplications over
    /// G1:
    /// - it prodceeds with `num_var` number of rounds,
    /// - at round i, we compute an MSM for `2^{num_var - i + 1}` number of G2
    ///   elements.
    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomial: &Self::Polynomial,
        point: &Self::Point,
        transcript: &mut IOPTranscript<C::ScalarField>,
    ) -> Result<(Self::Proof, C::ScalarField), PCSError> {
        let timer = start_timer!(|| "MultilinearHyraxPCS::open");
        let mut random_tape = RandomTape::new(b"proof");
        let eval = evaluate_opt(&polynomial, point);
        let proof_eval_vars_at_r = MultilinearHyraxProof::prove(
            &polynomial,
            None,
            &point,
            &eval,
            None,
            &prover_param.borrow(),
            transcript,
            &mut random_tape,
        )?;
        end_timer!(timer);
        Ok((proof_eval_vars_at_r, eval))
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    ///
    /// This function takes
    /// - num_var number of pairing product.
    /// - num_var number of MSM
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &Self::Point,
        value: &C::ScalarField,
        proof: &Self::Proof,
        transcript: &mut IOPTranscript<C::ScalarField>,
    ) -> Result<bool, PCSError> {
        Ok(proof.verify(verifier_param, transcript, point, value, commitment)?)
    }

    fn multi_open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomials: &[&[Self::Polynomial]],
        points: &[Self::Point],
        evals: &[&[C::ScalarField]],
        transcript: &mut IOPTranscript<C::ScalarField>,
    ) -> Result<Self::BatchProof, PCSError> {
        let timer = start_timer!(|| "MultilinearHyraxPCS::multi_open");
        let prover_param = prover_param.borrow();

        for point in points.iter() {
            // check whether point[i] == point[0]
            if *point != points[0] {
                return Err(PCSError::NotImplemented(String::from(
                    "Not implemented for different points",
                )));
            }
        }
        let proof = multi_open_single_point_internal(
            prover_param,
            polynomials,
            &points[0],
            evals,
            transcript,
        )?;

        end_timer!(timer);
        Ok(proof)
    }

    fn batch_verify(
        verifier_param: &Self::VerifierParam,
        commitments: &[Self::Commitment],
        points: &[Self::Point],
        batch_proof: &Self::BatchProof,
        transcript: &mut IOPTranscript<C::ScalarField>,
    ) -> Result<bool, PCSError> {
        for point in points.iter() {
            // check whether point[i] == point[0]
            if *point != points[0] {
                return Err(PCSError::NotImplemented(String::from(
                    "Not implemented for different points",
                )));
            }
        }
        batch_verify_single_point_internal(
            verifier_param,
            commitments,
            &points[0],
            batch_proof,
            transcript,
        )
    }
}

impl<C: AffineCurve> MultilinearHyraxPCS<C> {
    #[cfg(feature = "parallel")]
    fn commit_inner(
        poly: &DenseMultilinearExtension<C::ScalarField>,
        blinds: &[C::ScalarField],
        gens: &MultiCommitGens<C>,
    ) -> HyraxCommitment<C> {
        use rayon::prelude::{IntoParallelIterator, ParallelIterator};

        let (L_size, R_size) = (blinds.len(), gens.n);
        let comm = C::Projective::batch_normalization_into_affine(
            (0..L_size)
                .into_par_iter()
                .map(|i| {
                    commit_array(
                        &poly.evaluations[R_size * i..R_size * (i + 1)],
                        &blinds[i],
                        gens,
                    )
                })
                .collect::<Vec<_>>()
                .as_slice(),
        );
        HyraxCommitment { commitment: comm }
    }

    #[cfg(not(feature = "parallel"))]
    fn commit_inner(
        poly: &DenseMultilinearExtension<C::ScalarField>,
        blinds: &[C::ScalarField],
        gens: &MultiCommitGens<C>,
    ) -> HyraxCommitment<C> {
        let (left_num_vars, right_num_vars) =
            EqPolynomial::<C::ScalarField>::compute_factored_lens(poly.num_vars);
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();
        assert_eq!(L_size * R_size, poly.evaluations.len());
        let com = C::Projective::batch_normalization_into_affine(
            (0..L_size)
                .map(|i| {
                    commit_array(
                        &poly.evaluations[R_size * i..R_size * (i + 1)],
                        &blinds[i],
                        gens,
                    )
                })
                .collect::<Vec<_>>()
                .as_slice(),
        );
        HyraxCommitment { commitment: com }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::G1Affine as C;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::{rand::RngCore, test_rng, vec::Vec, UniformRand};

    type F = <C as AffineCurve>::ScalarField;

    fn test_single_helper<R: RngCore + CryptoRng>(
        params: &DotProductProofGens<C>,
        poly: &Arc<DenseMultilinearExtension<F>>,
        rng: &mut R,
    ) -> Result<(), PCSError> {
        let nv = poly.num_vars();
        assert_ne!(nv, 0);
        let (ck, vk) = MultilinearHyraxPCS::trim(params, None, Some(nv))?;
        let point: Vec<_> = (0..nv).map(|_| F::rand(rng)).collect();
        let com = MultilinearHyraxPCS::commit(&ck, poly)?;
        let mut transcript = IOPTranscript::<F>::new(b"test");
        let (proof, value) = MultilinearHyraxPCS::open(&ck, poly, &point, &mut transcript)?;
        let mut transcript = IOPTranscript::<F>::new(b"test");
        assert!(MultilinearHyraxPCS::verify(
            &vk,
            &com,
            &point,
            &value,
            &proof,
            &mut transcript,
        )?);

        let mut transcript = IOPTranscript::<F>::new(b"test");
        let value = F::rand(rng);
        assert!(!MultilinearHyraxPCS::verify(
            &vk,
            &com,
            &point,
            &value,
            &proof,
            &mut transcript
        )?);

        Ok(())
    }

    #[test]
    fn test_single_commit() -> Result<(), PCSError> {
        let mut rng = test_rng();

        let srs = MultilinearHyraxPCS::<C>::gen_srs_for_testing(&mut rng, 10)?;

        // normal polynomials
        let poly1 = Arc::new(DenseMultilinearExtension::rand(8, &mut rng));
        test_single_helper(&srs, &poly1, &mut rng)?;

        // single-variate polynomials
        let poly2 = Arc::new(DenseMultilinearExtension::rand(1, &mut rng));
        test_single_helper(&srs, &poly2, &mut rng)?;

        Ok(())
    }

    #[test]
    fn setup_commit_verify_constant_polynomial() {
        let mut rng = test_rng();

        // normal polynomials
        assert!(MultilinearHyraxPCS::<C>::gen_srs_for_testing(&mut rng, 0).is_err());
    }
}
