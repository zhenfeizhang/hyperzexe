//! Main module for the Lookup Check protocol

use crate::{
    pcs::PolynomialCommitmentScheme,
    poly_iop::{
        errors::PolyIOPErrors,
        lookup_check::util::{
            compute_count_in_input, compute_denominators_and_entry_frac, prove_zero_check,
        },
        prelude::ZeroCheck,
        PolyIOP,
    },
    SumCheck,
};
use arithmetic::VirtualPolynomial;
use ark_ec::AffineCurve;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer, One, Zero};
use std::sync::Arc;
use transcript::IOPTranscript;

/// A lookup subclaim consists of
/// - the SubClaim from the Sumcheck
/// - Challenges beta
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LookupCheckSubClaim<F: PrimeField, ZC: ZeroCheck<F>, SC: SumCheck<F>> {
    /// the SubClaim from the ZeroCheck
    pub zero_check_sub_claim: ZC::ZeroCheckSubClaim,
    /// the SubClaim from the Sumcheck
    pub sumcheck_sub_claim: SC::SumCheckSubClaim,
    /// Challenges beta
    pub beta: F,
}

/// A lookup check proof consists of
/// - a zerocheck proof
/// - a count polynomial commitment
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LookupCheckProof<
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<C, Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>>,
    ZC: ZeroCheck<C::ScalarField>,
    SC: SumCheck<C::ScalarField>,
> {
    pub zero_check_proof: ZC::ZeroCheckProof,
    pub sumcheck_proof: SC::SumCheckProof,
    pub count_comm: PCS::Commitment,
    pub frac_comm: PCS::Commitment,
}

pub mod util;

/// A LookupCheck w.r.t. `(fxs, gx)`
/// proves that all elements in f1, ..., fk belongs to table gx, where
/// A Lookup Check IOP takes the following steps:
///
/// Inputs:
/// - fxs = (f1, ..., fk)
/// - gx
///
/// Prover steps:
/// 1. build MLE `count(x)` as the number of occurrence of `gx` in `fxs`
/// 2. push commitments of `count(x)`
/// 3. Build zero check proof for virtual commitment
///     Q(x): \sum_(f \in fs) \prod_{f' \neq f}(f'(x) + \beta)(g(x) + \beta) -
/// m(x)\prod_{f \in fxs}(f(x) + \beta)         return a zero check proof for
/// Q(x)
///
/// Verifier steps:
/// 1. verify zero check proof
/// 2. verify polynomial openings
pub trait LookupCheck<C, PCS>: ZeroCheck<C::ScalarField>
where
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<C, Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>>,
{
    type LookupCheckSubClaim;
    type LookupProof;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a LookupCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// LookupCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Inputs:
    /// - inputs: fxs = (f1, ..., fk)
    /// - table: gx
    /// Outputs:
    /// - a lookup check proof proving that input fxs are subsets of table gx
    /// - the count polynomial built during lookup check
    /// Cost: O(N)
    #[allow(clippy::type_complexity)]
    fn prove(
        pcs_param: &PCS::ProverParam,
        fxs: &[Self::MultilinearExtension],
        gx: &Self::MultilinearExtension,
        transcript: &mut IOPTranscript<C::ScalarField>,
    ) -> Result<
        (
            Self::LookupProof,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
        ),
        PolyIOPErrors,
    >;

    /// Verify that (f1, ..., fk) are subsets of gx
    fn verify(
        proof: &Self::LookupProof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::LookupCheckSubClaim, PolyIOPErrors>;
}

impl<C, PCS> LookupCheck<C, PCS> for PolyIOP<C::ScalarField>
where
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<C, Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>>,
{
    type LookupCheckSubClaim = LookupCheckSubClaim<C::ScalarField, Self, Self>;
    type LookupProof = LookupCheckProof<C, PCS, Self, Self>;

    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<C::ScalarField>::new(b"Initializing PermutationCheck transcript")
    }

    fn prove(
        pcs_param: &PCS::ProverParam,
        fxs: &[Self::MultilinearExtension],
        gx: &Self::MultilinearExtension,
        transcript: &mut IOPTranscript<C::ScalarField>,
    ) -> Result<
        (
            Self::LookupProof,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
        ),
        PolyIOPErrors,
    > {
        let start = start_timer!(|| "Lookup check prove");
        if fxs.is_empty() {
            return Err(PolyIOPErrors::InvalidParameters(
                "input is empty".to_string(),
            ));
        }

        let count_poly = compute_count_in_input(fxs, gx)?;
        let count_comm = PCS::commit(pcs_param, &count_poly)?;
        transcript.append_serializable_element(b"count(x)", &count_comm)?;

        let beta = transcript.get_and_append_challenge(b"lookup beta")?;
        let (denominators, frac_poly) =
            compute_denominators_and_entry_frac(&fxs, &gx, &count_poly, beta)?;
        let frac_comm = PCS::commit(pcs_param, &frac_poly)?;
        transcript.append_serializable_element(b"lookup frac(x)", &frac_comm)?;

        // invoke correctness check on the entry_frac_poly
        let (zero_check_proof, _) =
            prove_zero_check(denominators.as_slice(), &frac_poly, &count_poly, transcript)?;

        // invoke sumcheck on entry_frac_poly
        let sum = VirtualPolynomial::new_from_mle(&frac_poly, C::ScalarField::one());
        let sumcheck_proof = <PolyIOP<_> as SumCheck<C::ScalarField>>::prove(&sum, transcript)?;

        end_timer!(start);
        Ok((
            Self::LookupProof {
                zero_check_proof,
                sumcheck_proof,
                count_comm,
                frac_comm,
            },
            count_poly,
            frac_poly,
        ))
    }

    fn verify(
        proof: &Self::LookupProof,
        zc_aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::LookupCheckSubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "prod_check verify");

        // update transcript and generate challenge
        transcript.append_serializable_element(b"count(x)", &proof.count_comm)?;
        let beta = transcript.get_and_append_challenge(b"lookup beta")?;
        transcript.append_serializable_element(b"lookup frac(x)", &proof.frac_comm)?;
        let zero_check_sub_claim = <Self as ZeroCheck<C::ScalarField>>::verify(
            &proof.zero_check_proof,
            zc_aux_info,
            transcript,
        )?;

        let mut sc_aux_info = zc_aux_info.clone();
        sc_aux_info.max_degree = 1;
        let sumcheck_sub_claim = <Self as SumCheck<C::ScalarField>>::verify(
            C::ScalarField::zero(),
            &proof.sumcheck_proof,
            &sc_aux_info,
            transcript,
        )?;

        end_timer!(start);

        Ok(LookupCheckSubClaim {
            zero_check_sub_claim,
            sumcheck_sub_claim,
            beta,
        })
    }
}

#[cfg(test)]
mod test {
    use super::LookupCheck;
    use crate::{
        pcs::{prelude::MultilinearHyraxPCS, PolynomialCommitmentScheme},
        poly_iop::{errors::PolyIOPErrors, PolyIOP},
    };
    use arithmetic::VPAuxInfo;
    use ark_bls12_381::{Fr, G1Affine as G1};
    use ark_ec::AffineCurve;
    use ark_ff::PrimeField;
    use ark_poly::DenseMultilinearExtension;
    use ark_std::test_rng;
    use std::{marker::PhantomData, sync::Arc};

    fn check_count_poly<C: AffineCurve>(
        count_poly: &Arc<DenseMultilinearExtension<C::ScalarField>>,
        fs: &[Arc<DenseMultilinearExtension<C::ScalarField>>],
        g: &Arc<DenseMultilinearExtension<C::ScalarField>>,
    ) -> Result<(), PolyIOPErrors> {
        let mut evals = vec![];
        for (g_eval, count_eval) in g.iter().zip(count_poly.iter()) {
            let count = count_eval.into_repr();
            evals.extend(vec![
                g_eval.into_repr().as_ref()[0] as u64;
                count.as_ref()[0] as usize
            ]);
        }
        let mut fs_evals: Vec<_> = fs
            .iter()
            .flat_map(|f| {
                f.iter()
                    .map(|f_eval| f_eval.into_repr().as_ref()[0] as u64)
                    .collect::<Vec<_>>()
            })
            .collect();

        // sort evals and fs_evals
        evals.sort();
        fs_evals.sort();

        assert_eq!(evals, fs_evals);
        Ok(())
    }

    fn test_lookup_check_helper<C, PCS>(
        pcs_param: &PCS::ProverParam,
        fxs: &[Arc<DenseMultilinearExtension<C::ScalarField>>],
        gx: &Arc<DenseMultilinearExtension<C::ScalarField>>,
    ) -> Result<(), PolyIOPErrors>
    where
        C: AffineCurve,
        PCS: PolynomialCommitmentScheme<
            C,
            Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>,
        >,
    {
        let nv = fxs[0].num_vars;
        // what's AuxInfo used for?
        let poly_info = VPAuxInfo {
            max_degree: fxs.len() + 2,
            num_variables: nv,
            phantom: PhantomData::default(),
        };

        // prover
        let mut transcript = <PolyIOP<C::ScalarField> as LookupCheck<C, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;
        let (proof, count, _) = <PolyIOP<C::ScalarField> as LookupCheck<C, PCS>>::prove(
            pcs_param,
            fxs,
            gx,
            &mut transcript,
        )?;

        // verifier
        let mut transcript = <PolyIOP<C::ScalarField> as LookupCheck<C, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;
        let _ = <PolyIOP<C::ScalarField> as LookupCheck<C, PCS>>::verify(
            &proof,
            &poly_info,
            &mut transcript,
        )?;

        check_count_poly::<C>(&count, fxs, gx)?;

        Ok(())
    }

    // write a function in rust, with two parameters fxs as a 2-D array and gx as an
    // array, check whether all elements in fxs belong to gx
    fn test_lookup_check(
        inputs: &[&[usize]],
        table: &[usize],
        nv: usize,
    ) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        let srs = MultilinearHyraxPCS::gen_srs_for_testing(&mut rng, nv)?;
        let (pcs_param, _) = MultilinearHyraxPCS::<G1>::trim(&srs, None, Some(nv))?;

        let fxs: Vec<_> = inputs
            .into_iter()
            .map(|input| {
                let fx: Vec<_> = input.into_iter().map(|x| Fr::from(*x as u64)).collect();
                Arc::new(DenseMultilinearExtension::from_evaluations_vec(nv, fx))
            })
            .collect();
        let table: Vec<_> = table.into_iter().map(|x| Fr::from(*x as u64)).collect();
        let gx = Arc::new(DenseMultilinearExtension::from_evaluations_vec(nv, table));
        test_lookup_check_helper::<G1, MultilinearHyraxPCS<G1>>(&pcs_param, &fxs, &gx)
    }

    #[test]
    fn test_normal_set_table_single_input() -> Result<(), PolyIOPErrors> {
        test_lookup_check(&[&[2, 2]], &[1, 2], 1)
    }

    #[test]
    fn test_normal_set_table_multi_input() -> Result<(), PolyIOPErrors> {
        test_lookup_check(&[&[1, 1], &[1, 2]], &[1, 2], 1)
    }

    #[test]
    fn test_multi_set_table() -> Result<(), PolyIOPErrors> {
        test_lookup_check(&[&[1, 1, 1, 1], &[1, 2, 2, 2]], &[1, 1, 2, 2], 2)
    }

    #[test]
    fn test_err_case() -> Result<(), PolyIOPErrors> {
        assert!(test_lookup_check(&[&[1, 2]], &[2, 3], 1).is_err());
        Ok(())
    }
}
