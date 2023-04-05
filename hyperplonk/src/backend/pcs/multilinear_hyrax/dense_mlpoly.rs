#![allow(clippy::too_many_arguments)]
#![allow(non_snake_case)]
use std::cmp::min;

use halo2_curves::{
    group::{ff::Field, Curve},
    CurveAffine,
};

use super::{
    commitments::commit_element,
    math::Math,
    nizk::{DotProductProofGens, DotProductProofLog},
    random::RandomTape,
};

use crate::backend::{
    pcs::multilinear_hyrax::HyraxCommitment,
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::{sum, variable_base_msm},
        transcript::TranscriptWrite,
    },
    Error,
};

#[derive(Debug, Clone)]
pub struct PolyCommitmentGens<C: CurveAffine> {
    pub max_num_vars: usize,
    pub right_num_vars: usize,
    pub gens: DotProductProofGens<C>,
}

impl<C: CurveAffine> PolyCommitmentGens<C> {
    // Create from existing SRS.
    pub fn new_from_srs(
        srs: &DotProductProofGens<C>,
        supported_num_vars: usize,
        supported_num_batches: usize,
    ) -> Result<PolyCommitmentGens<C>, Error> {
        let max_num_vars = supported_num_vars + supported_num_batches.next_power_of_two().log_2();
        // Set the right_num_vars to guarantee that the commitment of a single
        // polynomial won't be shorter than the length of inner product
        // argument.
        let right_num_vars = min(max_num_vars + 1 >> 1, supported_num_vars);
        let gens = srs.trim(1 << right_num_vars)?;
        debug_assert!(gens.n == 1 << right_num_vars);
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

pub struct PolyCommitmentBlinds<C: CurveAffine> {
    pub blinds: Vec<C::Scalar>,
}

pub struct EqPolynomial<F: Field> {
    r: Vec<F>,
}

impl<F: Field> EqPolynomial<F> {
    pub fn new(r: Vec<F>) -> Self {
        EqPolynomial { r }
    }

    pub fn compute_factored_evals(&self, right_num_vars: usize) -> (Vec<F>, Vec<F>) {
        let ell = self.r.len();

        let L = MultilinearPolynomial::eq_xy(&self.r[right_num_vars..ell]).into_evals();
        let R = MultilinearPolynomial::eq_xy(&self.r[..right_num_vars]).into_evals();

        (L, R)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultilinearHyraxProof<C: CurveAffine> {
    proof: DotProductProofLog<C>,
    blind_eval: C::Scalar,
}

impl<C: CurveAffine> MultilinearHyraxProof<C> {
    pub fn protocol_name() -> &'static [u8] {
        b"polynomial evaluation proof"
    }

    pub fn prove(
        poly: &MultilinearPolynomial<C::Scalar>,
        blinds_opt: Option<&PolyCommitmentBlinds<C>>,
        r: &[C::Scalar], // point at which the polynomial is evaluated
        Zr: &C::Scalar,  // evaluation of \widetilde{Z}(r)
        blind_Zr_opt: Option<&C::Scalar>, // specifies a blind for Zr
        gens: &PolyCommitmentGens<C>,
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
        random_tape: &mut RandomTape<C::Scalar>,
    ) -> Result<Self, Error> {
        // assert vectors are of the right size
        assert_eq!(poly.num_vars(), r.len());

        let (left_num_vars, right_num_vars) = gens.compute_factored_lens(r.len());
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();

        let default_blinds = PolyCommitmentBlinds {
            blinds: vec![C::Scalar::zero(); L_size],
        };
        let blinds = blinds_opt.map_or(&default_blinds, |p| p);

        debug_assert_eq!(blinds.blinds.len(), L_size);

        let zero = C::Scalar::zero();
        let blind_Zr = blind_Zr_opt.map_or(&zero, |p| p);

        // compute the L and R vectors
        let eq = EqPolynomial::new(r.to_vec());
        let (L, R) = eq.compute_factored_evals(right_num_vars);
        debug_assert_eq!(L.len(), L_size);
        debug_assert_eq!(R.len(), R_size);

        // compute the vector underneath L*Z and the L*blinds
        // compute vector-matrix product between L and Z viewed as a matrix
        let LZ = poly.fix_last_variables(&r[right_num_vars..]).into_evals();
        let LZ_blind: C::Scalar = sum((0..L.len()).map(|i| blinds.blinds[i] * L[i]));

        // a dot product proof of size R_size
        let (proof, _C_LR, _C_Zr_prime) = DotProductProofLog::prove(
            &gens.gens,
            transcript,
            random_tape,
            &LZ,
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
        transcript: &mut impl TranscriptWrite<C, C::Scalar>,
        r: &[C::Scalar], // point at which the polynomial is evaluated
        eval: &C::Scalar,
        comm: &HyraxCommitment<C>,
    ) -> Result<bool, Error> {
        // compute L and R
        let eq = EqPolynomial::new(r.to_vec());
        let (_, right_num_vars) = gens.compute_factored_lens(r.len());
        let (L, R) = eq.compute_factored_evals(right_num_vars);

        // compute a weighted sum of commitments and L
        let C_LZ = variable_base_msm(&L, comm).to_affine();
        let C_Zr = commit_element(eval, &self.blind_eval, &gens.gens.gens_1).to_affine();
        self.proof
            .verify(R.len(), &gens.gens, transcript, &R, &C_LZ, &C_Zr)
    }
}

#[cfg(test)]
mod tests {
    use crate::backend::util::{arithmetic::sum, test::std_rng};

    use super::*;
    use halo2_curves::bn256::Fr;

    fn evaluate_with_LR(Z: &[Fr], r: &[Fr]) -> Fr {
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
            .map(|i| sum((0..m).map(|j| L[j] * Z[j * m + i])))
            .collect::<Vec<Fr>>();

        // compute dot product between LZ and R
        sum((0..LZ.len()).map(|i| LZ[i] * R[i]))
    }

    #[test]
    fn check_polynomial_evaluation() {
        // Z = [1, 2, 1, 4]
        let Z = vec![Fr::one(), Fr::from(2), Fr::from(1), Fr::from(4)];

        // r = [3,4]
        let r = vec![Fr::from(3), Fr::from(4)];

        let eval_with_LR = evaluate_with_LR(&Z, &r);
        let poly = MultilinearPolynomial::new(Z);

        let eval = poly.evaluate(&r);
        assert_eq!(eval, Fr::from(28));
        assert_eq!(eval_with_LR, eval);
    }

    pub fn compute_factored_chis_at_r(r: &[Fr]) -> (Vec<Fr>, Vec<Fr>) {
        let mut L: Vec<Fr> = Vec::new();
        let mut R: Vec<Fr> = Vec::new();

        let ell = r.len();
        assert!(ell % 2 == 0); // ensure ell is even
        let n = ell.pow2();
        let m = n.square_root();

        // compute row vector L
        for i in 0..m {
            let mut chi_i = Fr::one();
            for j in 0..ell / 2 {
                let bit_j = ((m * i) & (1 << (r.len() - j - 1))) > 0;
                if bit_j {
                    chi_i *= r[j];
                } else {
                    chi_i *= Fr::one() - r[j];
                }
            }
            L.push(chi_i);
        }

        // compute column vector R
        for i in 0..m {
            let mut chi_i = Fr::one();
            for j in ell / 2..ell {
                let bit_j = (i & (1 << (r.len() - j - 1))) > 0;
                if bit_j {
                    chi_i *= r[j];
                } else {
                    chi_i *= Fr::one() - r[j];
                }
            }
            R.push(chi_i);
        }
        (L, R)
    }

    pub fn compute_chis_at_r(r: &[Fr]) -> Vec<Fr> {
        let ell = r.len();
        let n = ell.pow2();
        let mut chis: Vec<Fr> = Vec::new();
        for i in 0..n {
            let mut chi_i = Fr::one();
            for j in 0..r.len() {
                let bit_j = (i & (1 << (r.len() - j - 1))) > 0;
                if bit_j {
                    chi_i *= r[j];
                } else {
                    chi_i *= Fr::one() - r[j];
                }
            }
            chis.push(chi_i);
        }
        chis
    }

    pub fn compute_outerproduct(L: Vec<Fr>, R: Vec<Fr>) -> Vec<Fr> {
        assert_eq!(L.len(), R.len());
        (0..L.len())
            .map(|i| (0..R.len()).map(|j| L[i] * R[j]).collect::<Vec<Fr>>())
            .collect::<Vec<Vec<Fr>>>()
            .into_iter()
            .flatten()
            .collect::<Vec<Fr>>()
    }

    #[test]
    fn check_memoized_chis() {
        let mut csprng = std_rng();

        let s = 10;
        let mut r: Vec<Fr> = Vec::new();
        for _i in 0..s {
            r.push(Fr::random(&mut csprng));
        }
        let chis = tests::compute_chis_at_r(&r);
        let r_rev = {
            let mut r_rev = r.clone();
            r_rev.reverse();
            r_rev
        };
        let chis_m = MultilinearPolynomial::eq_xy(&r_rev).into_evals();
        assert_eq!(chis, chis_m);
    }

    #[test]
    fn check_factored_chis() {
        let mut csprng = std_rng();

        let s = 10;
        let mut r: Vec<Fr> = Vec::new();
        for _i in 0..s {
            r.push(Fr::random(&mut csprng));
        }
        let chis = MultilinearPolynomial::eq_xy(&r).into_evals();
        let right_num_vars = (r.len() + 1) / 2;
        let (L, R) = EqPolynomial::new(r).compute_factored_evals(right_num_vars);
        let O = compute_outerproduct(L, R);
        assert_eq!(chis, O);
    }

    #[test]
    fn check_memoized_factored_chis() {
        let mut csprng = std_rng();

        let s = 10;
        let mut r: Vec<Fr> = Vec::new();
        for _i in 0..s {
            r.push(Fr::random(&mut csprng));
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
