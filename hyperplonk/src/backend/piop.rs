use halo2_curves::group::{
    ff::{Field, PrimeField},
    prime::PrimeCurveAffine,
};
use itertools::Itertools;
use rand::RngCore;
use std::{fmt::Debug, hash::Hash, iter, marker::PhantomData, sync::Arc};

use crate::{
    backend::{
        piop::{
            preprocess::{compose, generate_preprocessed_poly_group},
            prover::*,
            verifier::verify_sum_check,
        },
        util::{end_timer, start_timer},
    },
    Error,
};

use self::{lookup_check::LookupCheck, perm_check::PermutationCheck, prover::instances_polys};
use super::{
    pcs::{prelude::HyraxBatchProof, PolynomialCommitmentScheme},
    poly::multilinear::MultilinearPolynomial,
    util::{expression::Expression, transcript::TranscriptWrite},
    PlonkishBackend, PlonkishCircuit, PlonkishCircuitInfo,
};
use itertools::max;

mod lookup_check;
mod perm_check;
mod preprocess;
mod prover;
mod sum_check;
mod verifier;

#[cfg(any(test, feature = "benchmark"))]
pub mod util;

#[derive(Clone, Debug)]
pub struct HyperPlonk<Pcs>(PhantomData<Pcs>);

#[derive(Clone, Debug)]
pub struct HyperPlonkProof<F, Pcs>
where
    F: Field,
    Pcs: PolynomialCommitmentScheme<F>,
{
    witness_comms: Vec<Pcs::Commitment>,
    lookup_count_comm: Pcs::Commitment,
    lookup_perm_intermediate_comm: Pcs::Commitment,
    constraint_sum_check_msgs: Vec<Vec<F>>,
    intermediate_evals: Vec<F>,
    opening_sum_check_msgs: Vec<Vec<F>>,
    final_evals: Vec<F>,
    batch_proof: Pcs::BatchProof,
}

#[derive(Clone, Debug)]
pub struct HyperPlonkProverParam<F, Pcs>
where
    F: Field,
    Pcs: PolynomialCommitmentScheme<F>,
{
    pcs: Pcs::ProverParam,
    num_instances: Vec<usize>,
    num_witness_polys: Vec<usize>,
    num_challenges: Vec<usize>,
    lookups: Vec<Vec<(Expression<F>, Expression<F>)>>,
    max_degree: usize,
    num_vars: usize,
    constraint_expression: Expression<F>,
    preprocess_polys: Vec<Arc<MultilinearPolynomial<F>>>,
    permutation_polys: Vec<(usize, Arc<MultilinearPolynomial<F>>)>,
    opening_expression: Expression<F>,
    segment_groups: Vec<Vec<(usize, usize)>>,
}

#[derive(Clone, Debug)]
pub struct HyperPlonkVerifierParam<F, Pcs>
where
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F>,
{
    pub pcs: Pcs::VerifierParam,
    pub num_instances: Vec<usize>,
    pub num_witness_polys: Vec<usize>,
    pub num_challenges: Vec<usize>,
    pub num_vars: usize,
    pub num_preprocess: usize,
    pub constraint_expression: Expression<F>,
    pub prep_perm_comm: Pcs::Commitment,
    pub open_expression: Expression<F>,
    segment_groups: Vec<Vec<(usize, usize)>>,
}

impl<F, C, Pcs> PlonkishBackend<F, Pcs> for HyperPlonk<Pcs>
where
    F: PrimeField + Ord + Hash,
    C: PrimeCurveAffine<Scalar = F>,
    Pcs: PolynomialCommitmentScheme<
        F,
        Curve = C,
        Polynomial = Arc<MultilinearPolynomial<F>>,
        Point = Vec<F>,
        Commitment = Vec<C>,
    >,
{
    type ProverParam = HyperPlonkProverParam<F, Pcs>;
    type VerifierParam = HyperPlonkVerifierParam<F, Pcs>;
    type Proof = HyperPlonkProof<F, Pcs>;

    fn setup(size: usize, rng: impl RngCore) -> Result<Pcs::SRS, Error> {
        Pcs::setup(rng, size << 1)
    }

    fn preprocess(
        param: &Pcs::SRS,
        mut circuit_info: PlonkishCircuitInfo<F>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), Error> {
        assert!(circuit_info.is_well_formed());
        circuit_info.initialize_permutation_info(circuit_info.permutation_polys().as_slice());
        let (preprocess_polys, permutation_polys) = generate_preprocessed_poly_group(&circuit_info);

        let num_vars = circuit_info.k;
        let (pcs_pp, pcs_vp) = {
            let batches = vec![
                circuit_info.preprocess_polys.len() + permutation_polys.len(),
                max(circuit_info.num_witness_polys.iter()).unwrap().clone(),
                circuit_info.lookups.len(),
                circuit_info.lookups.len() + circuit_info.num_permutation_chunks() * 2,
            ];
            let max_num_batches = max(batches.iter()).unwrap().clone();

            Pcs::trim(param, None, Some(num_vars), max_num_batches)?
        };

        let prep_perm_comm = {
            let polys = preprocess_polys
                .iter()
                .chain(permutation_polys.iter())
                .cloned()
                .collect_vec();
            Pcs::multi_commit(&pcs_pp, polys.as_slice())?
        };

        let segment_groups = {
            let mut offset = circuit_info.witness_offset();
            let mut witness_offsets = Vec::new();
            for num_witness_poly in circuit_info.num_witness_polys.iter() {
                witness_offsets.push(vec![(offset, *num_witness_poly)]);
                offset += *num_witness_poly;
            }
            let mut res = vec![vec![
                (
                    circuit_info.preprocess_offset(),
                    circuit_info.preprocess_polys.len(),
                ),
                (
                    circuit_info.permutation_offset(),
                    circuit_info.num_permutation_polys,
                ),
                (
                    circuit_info.permutation_offset(),
                    circuit_info.num_permutation_polys,
                ),
            ]];
            res.extend(witness_offsets);
            res.extend(vec![
                vec![(
                    circuit_info.lookup_count_offset(),
                    circuit_info.lookups.len(),
                )],
                vec![
                    (circuit_info.lookup_h_offset(), circuit_info.lookups.len()),
                    (
                        circuit_info.permutation_frac_offset(),
                        circuit_info.num_permutation_chunks(),
                    ),
                    (
                        circuit_info.permutation_prod_offset(),
                        circuit_info.num_permutation_chunks(),
                    ),
                ],
            ]);
            res
        };

        // Compose `VirtualPolynomialInfo`
        let (max_degree, constraint_expression, open_expression) = compose(&circuit_info);
        let vp = HyperPlonkVerifierParam {
            pcs: pcs_vp,
            num_instances: circuit_info.num_instances.clone(),
            num_witness_polys: circuit_info.num_witness_polys.clone(),
            num_challenges: circuit_info.num_challenges.clone(),
            num_vars,
            num_preprocess: preprocess_polys.len(),
            constraint_expression: constraint_expression.clone(),
            open_expression: open_expression.clone(),
            prep_perm_comm: prep_perm_comm.clone(),
            segment_groups: segment_groups.clone(),
        };
        let pp = HyperPlonkProverParam {
            pcs: pcs_pp,
            num_instances: circuit_info.num_instances.clone(),
            num_witness_polys: circuit_info.num_witness_polys.clone(),
            num_challenges: circuit_info.num_challenges.clone(),
            lookups: circuit_info.lookups.clone(),
            max_degree,
            num_vars,
            constraint_expression,
            opening_expression: open_expression,
            preprocess_polys,
            permutation_polys: circuit_info
                .permutation_polys()
                .into_iter()
                .zip(permutation_polys)
                .collect(),
            segment_groups,
        };
        Ok((pp, vp))
    }

    fn prove(
        pp: &Self::ProverParam,
        instances: &[&[F]],
        circuit: &impl PlonkishCircuit<F>,
        transcript: &mut impl TranscriptWrite<Pcs::Curve, F>,
        _: impl RngCore,
    ) -> Result<HyperPlonkProof<F, Pcs>, Error> {
        for (num_instances, instances) in pp.num_instances.iter().zip_eq(instances) {
            assert_eq!(instances.len(), *num_instances);
            for instance in instances.iter() {
                transcript.common_field_element(instance)?;
            }
        }
        let instances_polys = instances_polys(pp.num_vars, instances.iter().cloned());

        // Round 0..n
        let timer = start_timer(|| {
            format!(
                "witness_collector_for_{}_rounds",
                pp.num_witness_polys.len()
            )
        });
        let mut witness_comms = Vec::new();
        let mut witness_polys = Vec::new();
        let mut challenges = Vec::with_capacity(pp.num_challenges.iter().sum::<usize>() + 4);
        for (round, (num_witness_polys, num_challenges)) in pp
            .num_witness_polys
            .iter()
            .zip_eq(pp.num_challenges.iter())
            .enumerate()
        {
            let polys = circuit
                .synthesize(round, &challenges)?
                .into_iter()
                .map(|p| Arc::new(MultilinearPolynomial::new(p)))
                .collect_vec();
            assert_eq!(polys.len(), *num_witness_polys);

            witness_comms.push(Pcs::multi_commit_and_send(&pp.pcs, &polys, transcript)?);
            witness_polys.extend(polys);
            challenges.extend(transcript.squeeze_challenges(*num_challenges));
        }
        end_timer(timer);

        let permutation_polys = pp
            .permutation_polys
            .iter()
            .map(|(_, p)| p)
            .cloned()
            .collect_vec();

        let polys = iter::empty()
            .chain(instances_polys.iter())
            .chain(pp.preprocess_polys.iter())
            .chain(witness_polys.iter())
            .cloned()
            .collect_vec();

        // Round n

        let beta = transcript.squeeze_challenge();

        let timer = start_timer(|| format!("lookup_count_polys-{}", pp.lookups.len()));
        let mut lookup_check = LookupCheck::default();
        let lookup_slices = pp.lookups.iter().map(|l| l.as_slice()).collect_vec();
        let lookup_count_polys = lookup_check.commit_counts(
            lookup_slices.as_slice(),
            polys.as_slice(),
            challenges.as_slice(),
            &beta,
        )?;
        drop(lookup_slices);
        let lookup_count_comm =
            Pcs::multi_commit_and_send(&pp.pcs, lookup_count_polys.as_slice(), transcript)?;
        end_timer(timer);

        // Round n+1
        let gamma = transcript.squeeze_challenge();

        let timer = start_timer(|| format!("lookup_frac_polys-{}", pp.lookups.len()));
        let lookup_h_polys = lookup_check.commit_hs(lookup_count_polys.as_slice(), &gamma)?;

        end_timer(timer);

        let timer = start_timer(|| format!("permutation_z_polys-{}", pp.permutation_polys.len()));
        let (permutation_frac_polys, permutation_prod_polys, p1_p2_polys) =
            PermutationCheck::commit(pp.max_degree, &pp.permutation_polys, &polys, &beta, &gamma)?;

        let lookup_perm_intermediate_comm = {
            let to_commit_polys = lookup_h_polys
                .iter()
                .chain(permutation_frac_polys.iter())
                .chain(permutation_prod_polys.iter())
                .cloned()
                .collect_vec();
            Pcs::multi_commit_and_send(&pp.pcs, to_commit_polys.as_slice(), transcript)?
        };
        end_timer(timer);

        // Round n+2

        let alpha = transcript.squeeze_challenge();
        let y = transcript.squeeze_challenges(pp.num_vars);

        let p1_p2_polys = p1_p2_polys
            .into_iter()
            .flat_map(|(p1, p2)| vec![p1, p2])
            .collect_vec();

        // Compute all offsets to help reorder the polynomials to be consistent with the
        // commitments
        let preprocess_offset = instances_polys.len();
        let mut polys = polys;
        polys.extend(permutation_polys);
        polys.extend(lookup_count_polys);
        polys.extend(lookup_h_polys);
        polys.extend(permutation_frac_polys);
        polys.extend(permutation_prod_polys);

        let p1_p2_offset = polys.len();
        polys.extend(p1_p2_polys);

        challenges.extend([beta, gamma, alpha]);

        println!("Round n+2: sum_check-for-constraints");
        let timer = start_timer(|| format!("sum_check-for-constraints"));
        let (constraint_sum_check_msgs, points, intermediate_evals) = prove_sum_check(
            pp.num_instances.len(),
            &pp.constraint_expression,
            polys.as_slice(),
            challenges.as_slice(),
            y,
            F::zero(),
            transcript,
        )?;
        end_timer(timer);

        transcript.write_field_elements(&intermediate_evals)?;

        // Round n + 2 + num_vars
        let eta = transcript.squeeze_challenge();
        challenges.push(eta);

        println!("Round n+2+num_vars: sum_check-for-pc-openings");
        let timer = start_timer(|| format!("sum_check-for-pc-openings"));
        let random_combined_eval = intermediate_evals[preprocess_offset..]
            .iter()
            .fold(F::zero(), |acc, e| acc * eta + e);

        let (opening_sum_check_msgs, points, final_evals) = prove_sum_check(
            pp.num_instances.len(),
            &pp.opening_expression,
            &polys[0..p1_p2_offset],
            challenges.as_slice(),
            points[0].clone(),
            random_combined_eval,
            transcript,
        )?;
        end_timer(timer);
        transcript.write_field_elements(&final_evals[preprocess_offset..])?;

        // PCS open: reorder the polynomials and evaluations to match the order of the
        // commitments.
        let timer = start_timer(|| format!("pcs_batch_open"));

        let poly_groups = reorder_into_groups(polys, &pp.segment_groups);
        let eval_groups = reorder_into_groups(final_evals.clone(), &pp.segment_groups);

        let poly_group_slices = poly_groups
            .iter()
            .map(|polys| polys.as_slice())
            .collect_vec();
        let eval_group_slices = eval_groups
            .iter()
            .map(|evals| evals.as_slice())
            .collect_vec();
        let batch_proof = Pcs::multi_open(
            &pp.pcs,
            poly_group_slices.as_slice(),
            &points,
            eval_group_slices.as_slice(),
            transcript,
        )?;

        end_timer(timer);

        Ok(HyperPlonkProof {
            witness_comms,
            lookup_count_comm,
            lookup_perm_intermediate_comm,
            constraint_sum_check_msgs,
            intermediate_evals,
            opening_sum_check_msgs,
            final_evals,
            batch_proof,
        })
    }

    fn verify(
        vp: &Self::VerifierParam,
        instances: &[&[F]],
        proof: &Self::Proof,
        transcript: &mut impl TranscriptWrite<Pcs::Curve, F>,
        _: impl RngCore,
    ) -> Result<(), Error> {
        for (num_instances, instances) in vp.num_instances.iter().zip_eq(instances) {
            assert_eq!(instances.len(), *num_instances);
            for instance in instances.iter() {
                transcript.common_field_element(instance)?;
            }
        }

        // Round 0..n
        let mut challenges = Vec::with_capacity(vp.num_challenges.iter().sum::<usize>() + 4);
        for (witness_comm, num_challenges) in
            proof.witness_comms.iter().zip_eq(vp.num_challenges.iter())
        {
            transcript.write_commitments(witness_comm)?;
            challenges.extend(transcript.squeeze_challenges(*num_challenges));
        }

        // Round n
        let beta = transcript.squeeze_challenge();
        transcript.write_commitments(proof.lookup_count_comm.iter())?;

        // Round n+1
        let gamma = transcript.squeeze_challenge();
        transcript.write_commitments(proof.lookup_perm_intermediate_comm.iter())?;

        // Round n+2
        let alpha = transcript.squeeze_challenge();
        let y = transcript.squeeze_challenges(vp.num_vars);

        challenges.extend([beta, gamma, alpha]);
        let msg_slices = proof
            .constraint_sum_check_msgs
            .iter()
            .map(|v| v.as_slice())
            .collect_vec();
        let points = verify_sum_check(
            vp.num_vars,
            &vp.constraint_expression,
            instances,
            msg_slices.as_slice(),
            &challenges,
            &y,
            F::zero(),
            proof.intermediate_evals.as_slice(),
            transcript,
        )?;

        transcript.write_field_elements(&proof.intermediate_evals)?;

        // Round n + 2 + num_vars
        let eta = transcript.squeeze_challenge();
        challenges.push(eta);

        let msg_slices = proof
            .opening_sum_check_msgs
            .iter()
            .map(|v| v.as_slice())
            .collect_vec();
        let random_combined_eval = proof.intermediate_evals[instances.len()..]
            .iter()
            .fold(F::zero(), |acc, e| acc * eta + e);
        let points = verify_sum_check(
            vp.num_vars,
            &vp.open_expression,
            instances,
            msg_slices.as_slice(),
            &challenges,
            points[0].as_slice(),
            random_combined_eval,
            &proof.final_evals,
            transcript,
        )?;
        transcript.write_field_elements(&proof.final_evals[instances.len()..])?;
        let eval_groups = reorder_into_groups(proof.final_evals.clone(), &vp.segment_groups);
        let eval_group_slices = eval_groups
            .iter()
            .map(|evals| evals.as_slice())
            .collect_vec();

        // PCS verify
        let comms = iter::empty()
            .chain(iter::once(&vp.prep_perm_comm))
            .chain(&proof.witness_comms)
            .chain(iter::once(&proof.lookup_count_comm))
            .chain(iter::once(&proof.lookup_perm_intermediate_comm))
            .collect_vec();
        Pcs::batch_verify(
            &vp.pcs,
            &comms,
            &points,
            &eval_group_slices,
            &proof.batch_proof,
            transcript,
        )?;

        Ok(())
    }
}

#[cfg(test)]
pub(crate) mod test {
    use crate::backend::{
        pcs::prelude::MultilinearHyraxPCS,
        piop::{
            util::{rand_plonk_circuit, rand_plonk_with_lookup_circuit},
            HyperPlonk, PolynomialCommitmentScheme,
        },
        poly::multilinear::MultilinearPolynomial,
        util::{
            end_timer, start_timer,
            transcript::{
                InMemoryTranscriptRead, InMemoryTranscriptWrite, Keccak256Transcript,
                TranscriptRead, TranscriptWrite,
            },
            Itertools,
        },
        PlonkishBackend, PlonkishCircuit, PlonkishCircuitInfo,
    };
    use halo2_curves::{
        bn256::G1Affine as G1,
        group::{ff::PrimeField, prime::PrimeCurveAffine},
    };
    use rand::rngs::OsRng;
    use std::{hash::Hash, ops::Range, sync::Arc};

    pub(crate) fn run_hyperplonk<'a, F, Pcs, C, T, Circuit>(
        num_vars_range: Range<usize>,
        circuit_fn: impl Fn(usize) -> (PlonkishCircuitInfo<F>, Vec<Vec<F>>, Circuit),
    ) where
        F: PrimeField + Ord + Hash,
        C: PrimeCurveAffine<Scalar = F>,
        Pcs: PolynomialCommitmentScheme<
            F,
            Curve = C,
            Polynomial = Arc<MultilinearPolynomial<F>>,
            Point = Vec<F>,
            Commitment = Vec<C>,
        >,
        T: TranscriptRead<C, F>
            + TranscriptWrite<C, F>
            + InMemoryTranscriptRead
            + InMemoryTranscriptWrite,
        Circuit: PlonkishCircuit<F>,
    {
        for num_vars in num_vars_range {
            let (circuit_info, instances, circuit) = circuit_fn(num_vars);
            let instances = instances.iter().map(Vec::as_slice).collect_vec();

            let timer = start_timer(|| format!("setup-{num_vars}"));
            let param = HyperPlonk::<Pcs>::setup(1 << num_vars, OsRng).unwrap();
            end_timer(timer);

            let timer = start_timer(|| format!("preprocess-{num_vars}"));
            let (pp, vp) = HyperPlonk::<Pcs>::preprocess(&param, circuit_info).unwrap();
            end_timer(timer);

            let timer = start_timer(|| format!("prove-{num_vars}"));
            let proof = {
                let mut transcript = T::default();
                HyperPlonk::<Pcs>::prove(&pp, &instances, &circuit, &mut transcript, OsRng).unwrap()
            };
            end_timer(timer);

            let timer = start_timer(|| format!("verify-{num_vars}"));
            let result = {
                let mut transcript = T::default();
                HyperPlonk::<Pcs>::verify(&vp, &instances, &proof, &mut transcript, OsRng)
            };
            assert_eq!(result, Ok(()));
            end_timer(timer);
        }
    }

    macro_rules! tests {
        ($name:ident, $pcs:ty, $curve:ty) => {
            paste::paste! {
                #[test]
                fn [<$name _hyperplonk_plonk>]() {
                    run_hyperplonk::<_, $pcs, $curve, Keccak256Transcript<_>, _>(2..16, |num_vars| {
                        rand_plonk_circuit(num_vars, OsRng)
                    });
                }

                #[test]
                fn [<$name _hyperplonk_plonk_with_lookup>]() {
                    run_hyperplonk::<_, $pcs, $curve, Keccak256Transcript<_>, _>(2..16, |num_vars| {
                        rand_plonk_with_lookup_circuit(num_vars, OsRng)
                    });
                }
            }
        };
    }

    tests!(hyrax, MultilinearHyraxPCS<G1>, G1);
}
