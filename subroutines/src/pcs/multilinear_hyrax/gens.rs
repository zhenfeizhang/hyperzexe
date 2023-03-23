use super::nizk::DotProductProofGens;
use crate::PCSError;

pub struct PolyCommitmentGens {
    pub gens: DotProductProofGens,
}

impl PolyCommitmentGens {
    // the number of variables in the multilinear polynomial
    pub fn new(num_vars: usize, label: &'static [u8]) -> Result<PolyCommitmentGens, PCSError> {
        let (_left, right) = Self::compute_factored_lens(num_vars);
        let gens = DotProductProofGens::new(right.pow2(), label)?;
        Ok(PolyCommitmentGens { gens })
    }

    pub fn compute_factored_lens(ell: usize) -> (usize, usize) {
        (ell / 2, ell - ell / 2)
    }
}

pub struct PolyCommitmentBlinds {
    blinds: Vec<Scalar>,
}
