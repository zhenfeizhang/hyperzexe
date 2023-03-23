//! Main module for the HyperPlonk SNARK.

use std::sync::Arc;

use arithmetic::DenseMultilinearExtension;
use ark_ec::AffineCurve;
use errors::HyperPlonkErrors;
use subroutines::{pcs::prelude::PolynomialCommitmentScheme, poly_iop::prelude::PermutationCheck};
use witness::WitnessColumn;

mod custom_gate;
mod errors;
mod mock;
pub mod prelude;
mod selectors;
mod snark;
mod structs;
mod utils;
mod witness;

/// A trait for HyperPlonk SNARKs.
/// A HyperPlonk is derived from ZeroChecks and PermutationChecks.
pub trait HyperPlonkSNARK<C, PCS>: PermutationCheck<C, PCS>
where
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<C, Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>>,
{
    type Index;
    type ProvingKey;
    type VerifyingKey;
    type Proof;

    /// Generate the preprocessed polynomials output by the indexer.
    ///
    /// Inputs:
    /// - `index`: HyperPlonk index
    /// - `pcs_srs`: Polynomial commitment structured reference string
    /// Outputs:
    /// - The HyperPlonk proving key, which includes the preprocessed
    ///   polynomials.
    /// - The HyperPlonk verifying key, which includes the preprocessed
    ///   polynomial commitments
    fn preprocess(
        index: &Self::Index,
        pcs_srs: &PCS::SRS,
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), HyperPlonkErrors>;

    /// Generate HyperPlonk SNARK proof.
    ///
    /// Inputs:
    /// - `pk`: circuit proving key
    /// - `pub_input`: online public input
    /// - `witness`: witness assignment
    /// Outputs:
    /// - The HyperPlonk SNARK proof.
    fn prove(
        pk: &Self::ProvingKey,
        pub_input: &[C::ScalarField],
        witnesses: &[WitnessColumn<C::ScalarField>],
    ) -> Result<Self::Proof, HyperPlonkErrors>;

    /// Verify the HyperPlonk proof.
    ///
    /// Inputs:
    /// - `vk`: verifying key
    /// - `pub_input`: online public input
    /// - `proof`: HyperPlonk SNARK proof challenges
    /// Outputs:
    /// - Return a boolean on whether the verification is successful
    fn verify(
        vk: &Self::VerifyingKey,
        pub_input: &[C::ScalarField],
        proof: &Self::Proof,
    ) -> Result<bool, HyperPlonkErrors>;
}
