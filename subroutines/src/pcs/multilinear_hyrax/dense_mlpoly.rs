#![allow(clippy::too_many_arguments)]
#![allow(non_snake_case)]
use super::{
    commitments::commit_element,
    math::Math,
    nizk::{DotProductProofGens, DotProductProofLog},
    random::RandomTape,
};


use crate::{pcs::multilinear_hyrax::HyraxCommitment, PCSError};

use arithmetic::fix_last_variables;
use transcript::IOPTranscript;

use ark_ec::{msm::VariableBaseMSM, AffineCurve, ProjectiveCurve};
use ark_ff::{Field, PrimeField};
use ark_poly::DenseMultilinearExtension;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::Zero;

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct PolyCommitmentGens<C: AffineCurve> {
    pub max_num_vars: usize,
    pub right_num_vars: usize,
    pub gens: DotProductProofGens<C>,
}

impl<C: AffineCurve> PolyCommitmentGens<C> {
    // Create from the number of variables in the multilinear polynomial
    pub fn new(
        max_num_vars: usize,
        label: &'static [u8],
    ) -> Result<PolyCommitmentGens<C>, PCSError> {
        let right_num_vars = max_num_vars + 1 >> 1;
        let gens = DotProductProofGens::new(right_num_vars.pow2(), label)?;
        Ok(PolyCommitmentGens {
            max_num_vars,
            right_num_vars,
            gens,
        })
    }
    // Create from existing SRS.
    pub fn new_from_srs(
        srs: &DotProductProofGens<C>,
        max_num_vars: usize,
    ) -> Result<PolyCommitmentGens<C>, PCSError> {
        let right_num_vars = max_num_vars + 1 >> 1;
        let gens = srs.trim(1 << right_num_vars)?;
        Ok(PolyCommitmentGens {
            max_num_vars,
            right_num_vars,
            gens,
        })
    }
    pub fn compute_factored_lens(&self, num_vars: usize) -> (usize, usize) {
        (num_vars - self.right_num_vars, self.right_num_vars)
    }
}

pub struct PolyCommitmentBlinds<C: AffineCurve> {
    pub blinds: Vec<C::ScalarField>,
}

pub struct EqPolynomial<F: Field> {
    r: Vec<F>,
}

impl<F: Field> EqPolynomial<F> {
    pub fn new(r: Vec<F>) -> Self {
        EqPolynomial { r }
    }

    pub fn evals(&self) -> Vec<F> {
        let ell = self.r.len();

        let mut evals: Vec<F> = vec![F::one(); ell.pow2()];
        let mut size = 1;
        for j in 0..ell {
            // in each iteration, we double the size of chis
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                // copy each element from the prior iteration twice
                let scalar = evals[i / 2];
                evals[i] = scalar * self.r[ell - 1 - j];
                evals[i - 1] = scalar - evals[i];
            }
        }
        evals
    }

    pub fn compute_factored_evals(&self, right_num_vars: usize) -> (Vec<F>, Vec<F>) {
        let ell = self.r.len();

        let L = EqPolynomial::new(self.r[right_num_vars..ell].to_vec()).evals();
        let R = EqPolynomial::new(self.r[..right_num_vars].to_vec()).evals();

        (L, R)
    }
}

#[derive(Clone, Debug, CanonicalDeserialize, CanonicalSerialize, PartialEq, Eq)]
pub struct MultilinearHyraxProof<C: AffineCurve> {
    proof: DotProductProofLog<C>,
    blind_eval: C::ScalarField,
}

impl<C: AffineCurve> MultilinearHyraxProof<C> {
    pub fn protocol_name() -> &'static [u8] {
        b"polynomial evaluation proof"
    }

    pub fn prove(
        poly: &DenseMultilinearExtension<C::ScalarField>,
        blinds_opt: Option<&PolyCommitmentBlinds<C>>,
        r: &[C::ScalarField], // point at which the polynomial is evaluated
        Zr: &C::ScalarField,  // evaluation of \widetilde{Z}(r)
        blind_Zr_opt: Option<&C::ScalarField>, // specifies a blind for Zr
        gens: &PolyCommitmentGens<C>,
        transcript: &mut IOPTranscript<C::ScalarField>,
        random_tape: &mut RandomTape<C::ScalarField>,
    ) -> Result<Self, PCSError> {
        transcript.append_protocol_name(Self::protocol_name())?;

        // assert vectors are of the right size
        assert_eq!(poly.num_vars, r.len());

        let (left_num_vars, right_num_vars) = gens.compute_factored_lens(r.len());
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();

        let default_blinds = PolyCommitmentBlinds {
            blinds: vec![C::ScalarField::zero(); L_size],
        };
        let blinds = blinds_opt.map_or(&default_blinds, |p| p);

        assert_eq!(blinds.blinds.len(), L_size);

        let zero = C::ScalarField::zero();
        let blind_Zr = blind_Zr_opt.map_or(&zero, |p| p);

        // compute the L and R vectors
        let eq = EqPolynomial::new(r.to_vec());
        let (L, R) = eq.compute_factored_evals(right_num_vars);
        assert_eq!(L.len(), L_size);
        assert_eq!(R.len(), R_size);

        // compute the vector underneath L*Z and the L*blinds
        // compute vector-matrix product between L and Z viewed as a matrix
        let LZ = fix_last_variables(poly, &r[right_num_vars..]).evaluations;
        let LZ_blind: C::ScalarField = (0..L.len()).map(|i| blinds.blinds[i] * L[i]).sum();

        // a dot product proof of size R_size
        let (proof, _C_LR, _C_Zr_prime) = DotProductProofLog::prove(
            &gens.gens,
            transcript,
            random_tape,
            LZ.as_slice(),
            &LZ_blind,
            &R,
            Zr,
            blind_Zr,
            false, // we don't need blind property here.
        )?;

        Ok(Self {
            proof,
            blind_eval: *blind_Zr,
        })
    }

    pub fn verify(
        &self,
        gens: &PolyCommitmentGens<C>,
        transcript: &mut IOPTranscript<C::ScalarField>,
        r: &[C::ScalarField], // point at which the polynomial is evaluated
        eval: &C::ScalarField,
        comm: &HyraxCommitment<C>,
    ) -> Result<bool, PCSError> {
        transcript.append_protocol_name(Self::protocol_name())?;

        // compute L and R
        let eq = EqPolynomial::new(r.to_vec());
        let (_, right_num_vars) = gens.compute_factored_lens(r.len());
        let (L, R) = eq.compute_factored_evals(right_num_vars);

        // compute a weighted sum of commitments and L
        let L = L
            .into_iter()
            .map(|x| x.into_repr())
            .collect::<Vec<<C::ScalarField as PrimeField>::BigInt>>();
        let C_LZ = VariableBaseMSM::multi_scalar_mul(&comm.commitment, &L.as_slice()).into_affine();
        let C_Zr = commit_element(eval, &self.blind_eval, &gens.gens.gens_1).into_affine();
        self.proof
            .verify(R.len(), &gens.gens, transcript, &R, &C_LZ, &C_Zr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arithmetic::evaluate_opt;
    use ark_bls12_381::Fr;
    use ark_ff::{One, UniformRand};
    use ark_poly::DenseMultilinearExtension;

    use ark_std::test_rng;

    type Scalar = Fr;

    fn evaluate_with_LR(Z: &[Scalar], r: &[Scalar]) -> Scalar {
        let eq = EqPolynomial::new(r.to_vec());
        let right_num_vars = (r.len() + 1) / 2;
        let (L, R) = eq.compute_factored_evals(right_num_vars);

        let ell = r.len();
        // ensure ell is even
        assert!(ell % 2 == 0);
        // compute n = 2^\ell
        let n = ell.pow2();
        // compute m = sqrt(n) = 2^{\ell/2}
        let m = n.square_root();

        // compute vector-matrix product between L and Z viewed as a matrix
        let LZ = (0..m)
            .map(|i| (0..m).map(|j| L[j] * Z[j * m + i]).sum())
            .collect::<Vec<Scalar>>();

        // compute dot product between LZ and R
        (0..LZ.len()).map(|i| LZ[i] * R[i]).sum()
    }

    #[test]
    fn check_polynomial_evaluation() {
        // Z = [1, 2, 1, 4]
        let Z = vec![
            Scalar::one(),
            Scalar::from(2),
            Scalar::from(1),
            Scalar::from(4),
        ];

        // r = [3,4]
        let r = vec![Scalar::from(3), Scalar::from(4)];

        let eval_with_LR = evaluate_with_LR(&Z, &r);
        let poly = DenseMultilinearExtension::from_evaluations_slice(r.len(), &Z);

        let eval = evaluate_opt(&poly, &r);
        assert_eq!(eval, Scalar::from(28));
        assert_eq!(eval_with_LR, eval);
    }

    pub fn compute_factored_chis_at_r(r: &[Scalar]) -> (Vec<Scalar>, Vec<Scalar>) {
        let mut L: Vec<Scalar> = Vec::new();
        let mut R: Vec<Scalar> = Vec::new();

        let ell = r.len();
        assert!(ell % 2 == 0); // ensure ell is even
        let n = ell.pow2();
        let m = n.square_root();

        // compute row vector L
        for i in 0..m {
            let mut chi_i = Scalar::one();
            for j in 0..ell / 2 {
                let bit_j = ((m * i) & (1 << (r.len() - j - 1))) > 0;
                if bit_j {
                    chi_i *= r[j];
                } else {
                    chi_i *= Scalar::one() - r[j];
                }
            }
            L.push(chi_i);
        }

        // compute column vector R
        for i in 0..m {
            let mut chi_i = Scalar::one();
            for j in ell / 2..ell {
                let bit_j = (i & (1 << (r.len() - j - 1))) > 0;
                if bit_j {
                    chi_i *= r[j];
                } else {
                    chi_i *= Scalar::one() - r[j];
                }
            }
            R.push(chi_i);
        }
        (L, R)
    }

    pub fn compute_chis_at_r(r: &[Scalar]) -> Vec<Scalar> {
        let ell = r.len();
        let n = ell.pow2();
        let mut chis: Vec<Scalar> = Vec::new();
        for i in 0..n {
            let mut chi_i = Scalar::one();
            for j in 0..r.len() {
                let bit_j = (i & (1 << (r.len() - j - 1))) > 0;
                if bit_j {
                    chi_i *= r[j];
                } else {
                    chi_i *= Scalar::one() - r[j];
                }
            }
            chis.push(chi_i);
        }
        chis
    }

    pub fn compute_outerproduct(L: Vec<Scalar>, R: Vec<Scalar>) -> Vec<Scalar> {
        assert_eq!(L.len(), R.len());
        (0..L.len())
            .map(|i| (0..R.len()).map(|j| L[i] * R[j]).collect::<Vec<Scalar>>())
            .collect::<Vec<Vec<Scalar>>>()
            .into_iter()
            .flatten()
            .collect::<Vec<Scalar>>()
    }

    #[test]
    fn check_memoized_chis() {
        let mut csprng = test_rng();

        let s = 10;
        let mut r: Vec<Scalar> = Vec::new();
        for _i in 0..s {
            r.push(Scalar::rand(&mut csprng));
        }
        let chis = tests::compute_chis_at_r(&r);
        let r_rev = {
            let mut r_rev = r.clone();
            r_rev.reverse();
            r_rev
        };
        let chis_m = EqPolynomial::new(r_rev).evals();
        assert_eq!(chis, chis_m);
    }

    #[test]
    fn check_factored_chis() {
        let mut csprng = test_rng();

        let s = 10;
        let mut r: Vec<Scalar> = Vec::new();
        for _i in 0..s {
            r.push(Scalar::rand(&mut csprng));
        }
        let chis = EqPolynomial::new(r.clone()).evals();
        let right_num_vars = (r.len() + 1) / 2;
        let (L, R) = EqPolynomial::new(r).compute_factored_evals(right_num_vars);
        let O = compute_outerproduct(L, R);
        assert_eq!(chis, O);
    }

    #[test]
    fn check_memoized_factored_chis() {
        let mut csprng = test_rng();

        let s = 10;
        let mut r: Vec<Scalar> = Vec::new();
        for _i in 0..s {
            r.push(Scalar::rand(&mut csprng));
        }
        let (L, R) = tests::compute_factored_chis_at_r(&r);
        let r_rev = {
            let mut r_rev = r.clone();
            r_rev.reverse();
            r_rev
        };
        let eq = EqPolynomial::new(r_rev);
        let right_num_vars = (r.len() + 1) / 2;
        let (L2, R2) = eq.compute_factored_evals(right_num_vars);
        assert_eq!(L, L2);
        assert_eq!(R, R2);
    }
}
