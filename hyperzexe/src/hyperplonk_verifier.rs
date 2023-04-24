use halo2_ecc::fields::FieldChip;
use halo2_proofs::curves::CurveAffine;
use hyperplonk::backend::util::expression::{Expression, Query};
use serde::{Deserialize, Serialize};

use crate::{EcPoint, ScalarPoint};

pub(crate) use halo2_base::{halo2_proofs, halo2_proofs::halo2curves as halo2_curves};

mod pcs;
mod sumcheck;
mod system;
mod util;
mod verifier;

#[derive(Clone, Debug)]
pub enum Error {}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Protocol<C>
where
    C: CurveAffine,
{
    pub num_instance: Vec<usize>,
    pub num_witness: Vec<usize>,
    pub num_challenge: Vec<usize>,
    pub num_vars: usize,

    pub num_preprocess: usize,
    pub constraint_expression: Expression<ScalarPoint<C>>,
    pub opening_expression: Expression<ScalarPoint<C>>,
    pub prep_perm_comm: Vec<EcPoint<C>>,

    pub evaluations: Vec<Query>,
    // Minor customization
    pub transcript_initial_state: Option<ScalarPoint<C>>,
}
