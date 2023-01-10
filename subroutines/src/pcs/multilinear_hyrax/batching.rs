//! Sumcheck based batch opening and verify commitment.
// TODO: refactoring this code to somewhere else
// currently IOP depends on PCS because perm check requires commitment.
// The sumcheck based batch opening therefore cannot stay in the PCS repo --
// which creates a cyclic dependency.

use crate::pcs::PolynomialCommitmentScheme;
use ark_ec::AffineCurve;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BatchProof<C, PCS>
where
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<C>,
{
    /// f_i(point_i)
    pub f_i_eval_at_point_i: Vec<C::ScalarField>,
    pub f_i_proof_at_point_i: Vec<PCS::Proof>,
}
