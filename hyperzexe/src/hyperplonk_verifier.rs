use halo2_proofs::curves::CurveAffine;
use hyperplonk::backend::util::expression::{Expression, Query};

use crate::halo2_verifier::loader::Loader;

pub(crate) use halo2_base::halo2_proofs;
#[cfg(feature = "halo2-pse")]
pub(crate) use poseidon;
#[cfg(feature = "halo2-axiom")]
pub(crate) use poseidon_axiom as poseidon;

pub use poseidon::Spec as PoseidonSpec;
use serde::{Deserialize, Serialize};

mod pcs;
mod sumcheck;
mod system;
mod util;
mod verifier;

#[derive(Clone, Debug)]
pub enum Error {
    InvalidInstances,
    InvalidLinearization,
    InvalidQuery(Query),
    InvalidChallenge(usize),
    AssertionFailure(String),
    Transcript(std::io::ErrorKind, String),
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Protocol<C, L>
where
    C: CurveAffine,
    L: Loader<C>,
{
    pub num_instance: Vec<usize>,
    pub num_witness: Vec<usize>,
    pub num_challenge: Vec<usize>,
    pub num_vars: usize,
    pub max_degree: usize,
    pub num_all_polys: usize,
    // To reorder polynomials before calling PCS method
    pub segment_groups: Vec<Vec<usize>>,

    pub num_preprocess: usize,
    pub constraint_expression: Expression<L::LoadedScalar>,
    pub opening_expression: Expression<L::LoadedScalar>,
}
