//! Main module for multilinear Hyrax commitment scheme (based on implementation
//! in Spartan)
#![allow(non_snake_case)]

use std::{borrow::Borrow, marker::PhantomData, sync::Arc};

use halo2_curves::{group::ff::Field, CurveAffine};
use rand::RngCore;

use self::{
    batching::{batch_verify_single_point_internal, HyraxBatchProof},
    dense_mlpoly::{MultilinearHyraxProof, PolyCommitmentBlinds, PolyCommitmentGens},
    nizk::DotProductProofGens,
};
use crate::{
    backend::{
        pcs::{
            multilinear_hyrax::{
                batching::multi_open_single_point_internal,
                commitments::{commit_array, commit_element, MultiCommitGens},
                math::Math,
                random::RandomTape,
            },
            PolynomialCommitmentScheme,
        },
        poly::multilinear::MultilinearPolynomial,
        util::{end_timer, start_timer, transcript::TranscriptWrite},
    },
    Error,
};

pub(crate) mod batching;
mod commitments;
mod dense_mlpoly;
mod math;
mod nizk;
mod random;

/// Hyrax Polynomial Commitment Scheme on multilinear polynomials
#[derive(Clone, Debug, Default)]
pub struct MultilinearHyraxPCS<C: CurveAffine> {
    #[doc(hidden)]
    _curve: PhantomData<C>,
}

pub type HyraxSRS<C> = DotProductProofGens<C>;
pub type HyraxProverParam<C> = PolyCommitmentGens<C>;
pub type HyraxVerifierParam<C> = PolyCommitmentGens<C>;
pub type HyraxCommitment<C> = Vec<C>;

impl<C: CurveAffine> PolynomialCommitmentScheme<C::Scalar> for MultilinearHyraxPCS<C> {
    // Parameters
    type Curve = C;
    type ProverParam = HyraxProverParam<C>;
    type VerifierParam = HyraxVerifierParam<C>;
    /// Structured reference string
    type SRS = HyraxSRS<C>;
    // Polynomial and its associated types
    type Polynomial = Arc<MultilinearPolynomial<C::Scalar>>;
    type Point = Vec<C::Scalar>;
    // Commitments and proofs
    type Commitment = HyraxCommitment<C>;
    type Proof = MultilinearHyraxProof<C>;
    type BatchProof = HyraxBatchProof<C::Scalar, Self>;

    /// Build SRS for testing.
    ///
    /// - For univariate polynomials, `log_size` is the log of maximum degree.
    /// - For multilinear polynomials, `log_size` is the number of variables.
    ///
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn setup(_rng: impl RngCore, size: usize) -> Result<Self::SRS, Error> {
        if size == 0 {
            return Err(Error::InvalidPcsParam(
                "constant polynomial not supported".to_string(),
            ));
        }
        let log_size = size.log_2();
        let timer = start_timer(|| "MultilinearHyraxPCS::gen_srs_for_testing");
        let right_num_vars = log_size - 1;
        let res = DotProductProofGens::new(1 << right_num_vars, b"new srs")?;
        end_timer(timer);
        Ok(res)
    }

    /// Trim the universal parameters to specialize the public parameters.
    /// Input both `supported_log_degree` for univariate and
    /// `supported_num_vars` for multilinear.
    fn trim(
        srs: impl Borrow<Self::SRS>,
        supported_degree: Option<usize>,
        supported_num_vars: Option<usize>,
        supported_num_batches: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), Error> {
        let timer = start_timer(|| "MultilinearHyraxPCS::trim");
        assert!(supported_degree.is_none());
        let supported_num_vars = match supported_num_vars {
            Some(p) => p,
            None => {
                return Err(Error::InvalidPcsParam(
                    "multilinear should receive a num_var param".to_string(),
                ))
            },
        };
        let param = PolyCommitmentGens::new_from_srs(
            srs.borrow(),
            supported_num_vars,
            supported_num_batches,
        )?;
        end_timer(timer);

        Ok((param.clone(), param))
    }

    /// Generate a commitment for a polynomial.
    ///
    /// This function takes `2^num_vars` number of scalar multiplications over
    /// G1.
    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, Error> {
        let prover_param = prover_param.borrow();
        let timer = start_timer(|| "MultilinearHyraxPCS::commit");
        let n = poly.num_vars().pow2();
        let ell = poly.num_vars();

        let (left_num_vars, right_num_vars) = prover_param.compute_factored_lens(ell);
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();
        debug_assert_eq!(L_size * R_size, n);

        // We don't need hiding property here.
        let blinds: PolyCommitmentBlinds<C> = PolyCommitmentBlinds {
            blinds: vec![C::Scalar::zero(); L_size],
        };

        let res = Self::commit_inner(poly, &blinds.blinds, &prover_param.gens.gens_n);
        end_timer(timer);
        Ok(res)
    }

    /// Generate a commitment for polynomials
    fn multi_commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &[Self::Polynomial],
    ) -> Result<Self::Commitment, Error> {
        let prover_param = prover_param.borrow();
        let timer = start_timer(|| "MultilinearHyraxPCS::multi_commit");

        // Concatenate all polynomial evaluations and pad to the next power of 2.
        let mut poly_evals = Vec::new();
        for p in poly {
            poly_evals.extend(p.evals());
        }
        poly_evals.extend(vec![
            C::Scalar::zero();
            poly_evals.len().next_power_of_two() - poly_evals.len()
        ]);

        let num_vars = poly_evals.len().log_2();
        let poly = MultilinearPolynomial::new(poly_evals);

        let (left_num_vars, _) = prover_param.compute_factored_lens(num_vars);
        let L_size = left_num_vars.pow2();

        // We don't need hiding property here.
        let blinds: PolyCommitmentBlinds<C> = PolyCommitmentBlinds {
            blinds: vec![C::Scalar::zero(); L_size],
        };

        let res = Self::commit_inner(&poly, &blinds.blinds, &prover_param.gens.gens_n);
        end_timer(timer);
        Ok(res)
    }

    /// Generate a commitment for polynomials and send it to the verifier.
    fn multi_commit_and_send(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &[Self::Polynomial],
        transcript: &mut impl TranscriptWrite<Self::Curve, C::Scalar>,
    ) -> Result<Self::Commitment, Error> {
        if poly.len() == 0 {
            return Ok(vec![]);
        }
        let commitment = Self::multi_commit(prover_param, poly)?;
        transcript.write_commitments(&commitment)?;
        Ok(commitment)
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
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
    ) -> Result<(Self::Proof, C::Scalar), Error> {
        let timer = start_timer(|| "MultilinearHyraxPCS::open");
        let mut random_tape = RandomTape::new(b"proof");
        let eval = polynomial.evaluate(&point);
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
        end_timer(timer);
        Ok((proof_eval_vars_at_r, eval))
    }

    fn multi_open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomials: &[&[Self::Polynomial]],
        points: &[Self::Point],
        evals: &[&[C::Scalar]],
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
    ) -> Result<Self::BatchProof, Error> {
        let timer = start_timer(|| "MultilinearHyraxPCS::multi_open");
        let prover_param = prover_param.borrow();

        for point in points.iter() {
            // check whether point[i] == point[0]
            if *point != points[0] {
                return Err(Error::NotImplemented(String::from(
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

        end_timer(timer);
        Ok(proof)
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
        value: &C::Scalar,
        proof: &Self::Proof,
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
    ) -> Result<bool, Error> {
        Ok(proof.verify(verifier_param, transcript, point, value, commitment)?)
    }

    fn batch_verify(
        verifier_param: &Self::VerifierParam,
        commitments: &[&Self::Commitment],
        points: &[Self::Point],
        batch_proof: &Self::BatchProof,
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
    ) -> Result<bool, Error> {
        for point in points.iter() {
            // check whether point[i] == point[0]
            if *point != points[0] {
                return Err(Error::NotImplemented(String::from(
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

impl<C: CurveAffine> MultilinearHyraxPCS<C> {
    #[cfg(feature = "parallel")]
    fn commit_inner(
        poly: &MultilinearPolynomial<C::Scalar>,
        blinds: &[C::Scalar],
        gens: &MultiCommitGens<C>,
    ) -> HyraxCommitment<C> {
        use rayon::prelude::{IntoParallelIterator, ParallelIterator};

        let (L_size, R_size) = (blinds.len(), gens.n);
        (0..L_size)
            .into_par_iter()
            .map(|i| {
                commit_array(
                    &poly.evals()[R_size * i..R_size * (i + 1)],
                    &blinds[i],
                    gens,
                )
                .into()
            })
            .collect()
    }

    #[cfg(not(feature = "parallel"))]
    fn commit_inner(
        poly: &MultilinearPolynomial<C::Scalar>,
        blinds: &[C::Scalar],
        gens: &MultiCommitGens<C>,
    ) -> HyraxCommitment<C> {
        let (L_size, R_size) = (blinds.len(), gens.n);
        (0..L_size)
            .into_par_iter()
            .map(|i| {
                commit_array(
                    &poly.evals()[R_size * i..R_size * (i + 1)],
                    &blinds[i],
                    gens,
                )
                .into()
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::backend::util::{test::std_rng, transcript::Keccak256Transcript};

    use super::*;
    use halo2_curves::bn256::{Fr as F, G1Affine as C};
    use itertools::Itertools;

    fn test_single_helper(
        params: &DotProductProofGens<C>,
        poly: &Arc<MultilinearPolynomial<F>>,
        point: Vec<F>,
        eval: F,
    ) -> Result<(), Error> {
        let nv = poly.num_vars();
        assert_ne!(nv, 0);

        let (ck, vk) = MultilinearHyraxPCS::trim(params, None, Some(nv), 1)?;
        let com = MultilinearHyraxPCS::commit(&ck, poly)?;

        let mut transcript = Keccak256Transcript::<Vec<u8>>::default();
        let (proof, value) = MultilinearHyraxPCS::open(&ck, poly, &point, &mut transcript)?;
        assert!(eval == value);
        let mut transcript = Keccak256Transcript::<Vec<u8>>::default();
        assert!(MultilinearHyraxPCS::verify(
            &vk,
            &com,
            &point,
            &value,
            &proof,
            &mut transcript,
        )?);

        let mut transcript = Keccak256Transcript::<Vec<u8>>::default();
        let value = eval + F::one();
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
    fn test_single_one_variable_commit() -> Result<(), Error> {
        let mut rng = std_rng();

        let srs = MultilinearHyraxPCS::<C>::setup(&mut rng, 10)?;

        let poly2 = Arc::new(MultilinearPolynomial::rand(1, &mut rng));
        let point = vec![F::from(5u64)];
        let value = poly2[0] * (F::one() - point[0]) + poly2[1] * point[0];
        test_single_helper(&srs, &poly2, point, value)?;

        Ok(())
    }

    #[test]
    fn test_single_constant_commit() -> Result<(), Error> {
        let mut rng = std_rng();

        let srs = MultilinearHyraxPCS::<C>::setup(&mut rng, 10)?;

        let poly0 = {
            let ele = F::random(&mut rng);
            let evals = (0..16).map(|_| ele).collect::<Vec<_>>();
            Arc::new(MultilinearPolynomial::new(evals))
        };
        let point = (0..4).map(|_| F::random(&mut rng)).collect::<Vec<_>>();
        test_single_helper(&srs, &poly0, point, poly0[0])?;

        Ok(())
    }

    #[test]
    fn test_single_normal_commit() -> Result<(), Error> {
        let mut rng = std_rng();

        let srs = MultilinearHyraxPCS::<C>::setup(&mut rng, 10)?;

        let poly1 = Arc::new(MultilinearPolynomial::new(
            [1, 2, 3, 4]
                .iter()
                .map(|u| F::from(*u as u64))
                .collect_vec(),
        ));
        let point = vec![F::from(2u64), F::from(3u64)];
        let value = F::from(9u64);
        test_single_helper(&srs, &poly1, point, value)?;

        Ok(())
    }

    #[test]
    fn setup_commit_verify_constant_polynomial() {
        let mut rng = std_rng();

        // normal polynomials
        assert!(MultilinearHyraxPCS::<C>::setup(&mut rng, 0).is_err());
    }
}
