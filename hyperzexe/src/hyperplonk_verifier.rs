use halo2_proofs::curves::CurveAffine;
use hyperplonk::backend::util::expression::{Expression, Query};
use serde::{Deserialize, Serialize};

mod pcs;
mod sumcheck;
mod system;
mod util;
mod verifier;

#[derive(Clone, Debug)]
pub enum Error {}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Protocol<C, L = loader::native::NativeLoader>
where
    C: CurveAffine,
    L: loader::Loader<C>,
{
    pub num_instance: Vec<usize>,
    pub num_witness: Vec<usize>,
    pub num_challenge: Vec<usize>,
    pub num_vars: usize,

    pub num_preprocess: usize,
    pub constraint_expression: Expression<C::Scalar>,
    pub opening_expression: Expression<C::Scalar>,
    #[serde(bound(
        serialize = "L::LoadedScalar: Serialize",
        deserialize = "L::LoadedScalar: Deserialize<'de>"
    ))]
    pub prep_perm_comm: Vec<L::LoadedEcPoint>,

    pub evaluations: Vec<Query>,
    // Minor customization
    #[serde(bound(
        serialize = "L::LoadedScalar: Serialize",
        deserialize = "L::LoadedScalar: Deserialize<'de>"
    ))]
    pub transcript_initial_state: Option<L::LoadedScalar>,
}
