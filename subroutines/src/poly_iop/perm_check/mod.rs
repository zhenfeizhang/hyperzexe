//! Main module for the Permutation Check protocol

use self::util::computer_nums_and_denoms;
use crate::{
    pcs::PolynomialCommitmentScheme,
    poly_iop::{errors::PolyIOPErrors, prelude::ProductCheck, PolyIOP},
};
use ark_ec::AffineCurve;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use std::sync::Arc;
use transcript::IOPTranscript;

/// A permutation subclaim consists of
/// - the SubClaim from the ProductCheck
/// - Challenges beta and gamma
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PermutationCheckSubClaim<C, PCS, PC>
where
    C: AffineCurve,
    PC: ProductCheck<C, PCS>,
    PCS: PolynomialCommitmentScheme<C, Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>>,
{
    /// the SubClaim from the ProductCheck
    pub product_check_sub_claim: PC::ProductCheckSubClaim,
    /// Challenges beta and gamma
    pub challenges: (C::ScalarField, C::ScalarField),
}

pub mod util;

/// A PermutationCheck w.r.t. `(fs, gs, perms)`
/// proves that (g1, ..., gk) is a permutation of (f1, ..., fk) under
/// permutation `(p1, ..., pk)`
/// It is derived from ProductCheck.
///
/// A Permutation Check IOP takes the following steps:
///
/// Inputs:
/// - fs = (f1, ..., fk)
/// - gs = (g1, ..., gk)
/// - permutation oracles = (p1, ..., pk)
pub trait PermutationCheck<C, PCS>: ProductCheck<C, PCS>
where
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<C, Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>>,
{
    type PermutationCheckSubClaim;
    type PermutationProof;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Inputs:
    /// - fs = (f1, ..., fk)
    /// - gs = (g1, ..., gk)
    /// - permutation oracles = (p1, ..., pk)
    /// Outputs:
    /// - a permutation check proof proving that gs is a permutation of fs under
    ///   permutation
    /// - the product polynomial built during product check
    /// - the fractional polynomial built during product check
    ///
    /// Cost: O(N)
    #[allow(clippy::type_complexity)]
    fn prove(
        pcs_param: &PCS::ProverParam,
        fxs: &[Self::MultilinearExtension],
        gxs: &[Self::MultilinearExtension],
        perms: &[Self::MultilinearExtension],
        transcript: &mut IOPTranscript<C::ScalarField>,
    ) -> Result<
        (
            Self::PermutationProof,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
        ),
        PolyIOPErrors,
    >;

    /// Verify that (g1, ..., gk) is a permutation of
    /// (f1, ..., fk) over the permutation oracles (perm1, ..., permk)
    fn verify(
        proof: &Self::PermutationProof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationCheckSubClaim, PolyIOPErrors>;
}

impl<C, PCS> PermutationCheck<C, PCS> for PolyIOP<C::ScalarField>
where
    C: AffineCurve,
    PCS: PolynomialCommitmentScheme<C, Polynomial = Arc<DenseMultilinearExtension<C::ScalarField>>>,
{
    type PermutationCheckSubClaim = PermutationCheckSubClaim<C, PCS, Self>;
    type PermutationProof = Self::ProductCheckProof;

    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<C::ScalarField>::new(b"Initializing PermutationCheck transcript")
    }

    fn prove(
        pcs_param: &PCS::ProverParam,
        fxs: &[Self::MultilinearExtension],
        gxs: &[Self::MultilinearExtension],
        perms: &[Self::MultilinearExtension],
        transcript: &mut IOPTranscript<C::ScalarField>,
    ) -> Result<
        (
            Self::PermutationProof,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
        ),
        PolyIOPErrors,
    > {
        let start = start_timer!(|| "Permutation check prove");
        if fxs.is_empty() {
            return Err(PolyIOPErrors::InvalidParameters("fxs is empty".to_string()));
        }
        if (fxs.len() != gxs.len()) || (fxs.len() != perms.len()) {
            return Err(PolyIOPErrors::InvalidProof(format!(
                "fxs.len() = {}, gxs.len() = {}, perms.len() = {}",
                fxs.len(),
                gxs.len(),
                perms.len(),
            )));
        }

        let num_vars = fxs[0].num_vars;
        for ((fx, gx), perm) in fxs.iter().zip(gxs.iter()).zip(perms.iter()) {
            if (fx.num_vars != num_vars) || (gx.num_vars != num_vars) || (perm.num_vars != num_vars)
            {
                return Err(PolyIOPErrors::InvalidParameters(
                    "number of variables unmatched".to_string(),
                ));
            }
        }

        // generate challenge `beta` and `gamma` from current transcript
        let beta = transcript.get_and_append_challenge(b"beta")?;
        let gamma = transcript.get_and_append_challenge(b"gamma")?;
        let (numerators, denominators) = computer_nums_and_denoms(&beta, &gamma, fxs, gxs, perms)?;

        // invoke product check on numerator and denominator
        let (proof, prod_poly, frac_poly) = <Self as ProductCheck<C, PCS>>::prove(
            pcs_param,
            &numerators,
            &denominators,
            transcript,
        )?;

        end_timer!(start);
        Ok((proof, prod_poly, frac_poly))
    }

    fn verify(
        proof: &Self::PermutationProof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationCheckSubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "Permutation check verify");

        let beta = transcript.get_and_append_challenge(b"beta")?;
        let gamma = transcript.get_and_append_challenge(b"gamma")?;

        // invoke the zero check on the iop_proof
        let product_check_sub_claim =
            <Self as ProductCheck<C, PCS>>::verify(proof, aux_info, transcript)?;

        end_timer!(start);
        Ok(PermutationCheckSubClaim {
            product_check_sub_claim,
            challenges: (beta, gamma),
        })
    }
}

#[cfg(test)]
mod test {
    use super::PermutationCheck;
    use crate::{
        pcs::{prelude::MultilinearHyraxPCS, PolynomialCommitmentScheme},
        poly_iop::{errors::PolyIOPErrors, PolyIOP},
    };
    use arithmetic::{evaluate_opt, identity_permutation_mles, random_permutation_mles, VPAuxInfo};
    use ark_bls12_381::G1Affine as G1;
    use ark_ec::AffineCurve;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use std::{marker::PhantomData, sync::Arc};

    type Hyrax = MultilinearHyraxPCS<G1>;

    fn test_permutation_check_helper<C, PCS>(
        pcs_param: &PCS::ProverParam,
        fxs: &[Arc<DenseMultilinearExtension<C::ScalarField>>],
        gxs: &[Arc<DenseMultilinearExtension<C::ScalarField>>],
        perms: &[Arc<DenseMultilinearExtension<C::ScalarField>>],
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
            max_degree: fxs.len() + 1,
            num_variables: nv,
            phantom: PhantomData::default(),
        };

        // prover
        let mut transcript =
            <PolyIOP<C::ScalarField> as PermutationCheck<C, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;
        let (proof, prod_x, _frac_poly) =
            <PolyIOP<C::ScalarField> as PermutationCheck<C, PCS>>::prove(
                pcs_param,
                fxs,
                gxs,
                perms,
                &mut transcript,
            )?;

        // verifier
        let mut transcript =
            <PolyIOP<C::ScalarField> as PermutationCheck<C, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;
        let perm_check_sub_claim = <PolyIOP<C::ScalarField> as PermutationCheck<C, PCS>>::verify(
            &proof,
            &poly_info,
            &mut transcript,
        )?;

        // check product subclaim
        if evaluate_opt(
            &prod_x,
            &perm_check_sub_claim.product_check_sub_claim.final_query.0,
        ) != perm_check_sub_claim.product_check_sub_claim.final_query.1
        {
            return Err(PolyIOPErrors::InvalidVerifier("wrong subclaim".to_string()));
        };

        Ok(())
    }

    fn test_permutation_check(nv: usize) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        let srs = MultilinearHyraxPCS::gen_srs_for_testing(&mut rng, nv)?;
        let (pcs_param, _) = MultilinearHyraxPCS::<G1>::trim(&srs, None, Some(nv))?;
        let id_perms = identity_permutation_mles(nv, 2);

        {
            // good path: (w1, w2) is a permutation of (w1, w2) itself under the identify
            // map
            let ws = vec![
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
            ];
            // perms is the identity map
            test_permutation_check_helper::<G1, Hyrax>(&pcs_param, &ws, &ws, &id_perms)?;
        }

        {
            // good path: f = (w1, w2) is a permutation of g = (w2, w1) itself under a map
            let mut fs = vec![
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
            ];
            let gs = fs.clone();
            fs.reverse();
            // perms is the reverse identity map
            let mut perms = id_perms.clone();
            perms.reverse();
            test_permutation_check_helper::<G1, Hyrax>(&pcs_param, &fs, &gs, &perms)?;
        }

        {
            // bad path 1: w is a not permutation of w itself under a random map
            let ws = vec![
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
            ];
            // perms is a random map
            let perms = random_permutation_mles(nv, 2, &mut rng);

            assert!(
                test_permutation_check_helper::<G1, Hyrax>(&pcs_param, &ws, &ws, &perms).is_err()
            );
        }

        {
            // bad path 2: f is a not permutation of g under a identity map
            let fs = vec![
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
            ];
            let gs = vec![
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
                Arc::new(DenseMultilinearExtension::rand(nv, &mut rng)),
            ];
            // s_perm is the identity map

            assert!(
                test_permutation_check_helper::<G1, Hyrax>(&pcs_param, &fs, &gs, &id_perms)
                    .is_err()
            );
        }

        Ok(())
    }

    #[test]
    fn test_trivial_polynomial() -> Result<(), PolyIOPErrors> {
        test_permutation_check(1)
    }
    #[test]
    fn test_normal_polynomial() -> Result<(), PolyIOPErrors> {
        test_permutation_check(5)
    }

    #[test]
    fn zero_polynomial_should_error() -> Result<(), PolyIOPErrors> {
        assert!(test_permutation_check(0).is_err());
        Ok(())
    }
}
