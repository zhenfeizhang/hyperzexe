mod multilinear_hyrax;
mod structs;

pub mod prelude;

use crate::{backend::util::transcript::TranscriptWrite, Error};
use halo2_curves::{group::ff::Field, CurveAffine};
use rand::RngCore;
use std::{borrow::Borrow, fmt::Debug};

/// This trait defines APIs for polynomial commitment schemes.
/// Note that for our usage of PCS, we do not require the hiding property.
pub trait PolynomialCommitmentScheme<F: Field>: Clone + Debug {
    /// Curve type
    type Curve: CurveAffine;
    /// Prover parameters
    type ProverParam: Clone + Debug + Sync;
    /// Verifier parameters
    type VerifierParam: Clone + Debug;
    /// Structured reference string
    type SRS: Clone + Debug;
    /// Polynomial and its associated types
    type Polynomial: Clone + Debug + PartialEq + Eq;
    /// Polynomial input domain
    type Point: Clone + Ord + Debug + Sync + PartialEq + Eq;
    /// Commitments
    type Commitment: Clone + Debug + PartialEq + Eq + Default;
    /// Proofs
    type Proof: Clone + Debug + PartialEq + Eq;
    /// Batch proofs
    type BatchProof: Clone + Debug;

    /// Build SRS for testing.
    ///
    /// - For univariate polynomials, `supported_size` is the maximum degree.
    /// - For multilinear polynomials, `supported_size` is the number of
    ///   variables.
    ///
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn setup(rng: impl RngCore, log_size: usize) -> Result<Self::SRS, Error>;

    /// Trim the universal parameters to specialize the public parameters.
    /// Input both `supported_degree` for univariate and
    /// `supported_num_vars` for multilinear.
    /// ## Note on function signature
    /// Usually, data structure like SRS and ProverParam are huge and users
    /// might wish to keep them in heap using different kinds of smart pointers
    /// (instead of only in stack) therefore our `impl Borrow<_>` interface
    /// allows for passing in any pointer type, e.g.: `trim(srs: &Self::SRS,
    /// ..)` or `trim(srs: Box<Self::SRS>, ..)` or `trim(srs: Arc<Self::SRS>,
    /// ..)` etc.
    fn trim(
        srs: impl Borrow<Self::SRS>,
        supported_degree: Option<usize>,
        supported_num_vars: Option<usize>,
        supported_num_batches: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), Error>;

    /// Generate a commitment for a polynomial
    /// ## Note on function signature
    /// Usually, data structure like SRS and ProverParam are huge and users
    /// might wish to keep them in heap using different kinds of smart pointers
    /// (instead of only in stack) therefore our `impl Borrow<_>` interface
    /// allows for passing in any pointer type, e.g.: `commit(prover_param:
    /// &Self::ProverParam, ..)` or `commit(prover_param:
    /// Box<Self::ProverParam>, ..)` or `commit(prover_param:
    /// Arc<Self::ProverParam>, ..)` etc.
    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, Error>;

    /// Generate a commitment for polynomials
    fn multi_commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &[Self::Polynomial],
    ) -> Result<Self::Commitment, Error>;

    /// Generate a commitment for polynomials and send it to the transcript.
    fn multi_commit_and_send(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &[Self::Polynomial],
        transcript: &mut impl TranscriptWrite<Self::Curve, F>,
    ) -> Result<Self::Commitment, Error>;

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same.
    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomial: &Self::Polynomial,
        point: &Self::Point,
        transcript: &mut impl TranscriptWrite<Self::Curve, F>,
    ) -> Result<(Self::Proof, F), Error>;

    /// Input a list of multilinear extensions, and a same number of points, and
    /// a transcript, compute a multi-opening for all the polynomials.
    fn multi_open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomials: &[&[Self::Polynomial]],
        points: &[Self::Point],
        evals: &[&[F]],
        transcript: &mut impl TranscriptWrite<Self::Curve, F>,
    ) -> Result<Self::BatchProof, Error>;

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &Self::Point,
        value: &F,
        proof: &Self::Proof,
        transcript: &mut impl TranscriptWrite<Self::Curve, F>,
    ) -> Result<bool, Error>;

    /// Verifies that `value_i` is the evaluation at `x_i` of the polynomial
    /// `poly_i` committed inside `comm`.
    fn batch_verify(
        verifier_param: &Self::VerifierParam,
        commitments: &[&Self::Commitment],
        points: &[Self::Point],
        values: &[&[F]],
        batch_proof: &Self::BatchProof,
        transcript: &mut impl TranscriptWrite<Self::Curve, F>,
    ) -> Result<bool, Error>;
}

#[derive(Clone, Debug)]
pub struct Evaluation<F: Field> {
    poly: usize,
    point: usize,
    value: F,
}

impl<F: Field> Evaluation<F> {
    pub fn new(poly: usize, point: usize, value: F) -> Self {
        Self { poly, point, value }
    }

    pub fn poly(&self) -> usize {
        self.poly
    }

    pub fn point(&self) -> usize {
        self.point
    }

    pub fn value(&self) -> &F {
        &self.value
    }
}
