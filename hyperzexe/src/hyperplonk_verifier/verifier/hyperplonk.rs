use std::{marker::PhantomData, ops::Mul};

use halo2_proofs::curves::CurveAffine;

use super::HyperPlonkVerifier;
use crate::{
    halo2_verifier::{loader::Loader, util::transcript::TranscriptRead},
    hyperplonk_verifier::{
        pcs::{MultiOpenScheme, OpenScheme},
        sumcheck::{SumcheckRoundVerifier, SumcheckVerifier},
        Error, Protocol,
    },
};

pub struct HyperPlonk<C, L, SRC, SC, OS, MOS>(PhantomData<(C, L, SRC, SC, OS, MOS)>);

impl<C, L, SRC, SC, OS, MOS> HyperPlonkVerifier<C, L, SRC, SC, OS, MOS>
    for HyperPlonk<C, L, SRC, SC, OS, MOS>
where
    C: CurveAffine,
    L: Loader<C>,
    SRC: SumcheckRoundVerifier<C, L>,
    SC: SumcheckVerifier<C, L, SRC>,
    OS: OpenScheme<C, L>,
    MOS: MultiOpenScheme<C, L, OS>,
{
    type VerifyingKey = HyperPlonkVerifyingKey<C, L, OS, MOS>;
    type Proof = HyperPlonkProof<C, L, SRC, SC, OS, MOS>;
    type Output = ();

    fn read_proof<T>(
        vk: &Self::VerifyingKey,
        protocol: &Protocol<C, L>,
        instances: &[Vec<L::LoadedScalar>],
        transcript: &mut T,
    ) -> Result<Self::Proof, Error>
    where
        T: TranscriptRead<C, L>,
    {
        debug_assert_eq!(
            protocol.num_instance,
            instances
                .iter()
                .map(|instances| instances.len())
                .collect_vec(),
            "Invalid Instances"
        );

        for instances in instances.iter() {
            for instance in instances.iter() {
                transcript.common_scalar(instance).unwrap();
            }
        }

        let (witnesses, challenges) = {
            let (witnesses, challenges): (Vec<_>, Vec<_>) = protocol
                .num_witness
                .iter()
                .zip(protocol.num_challenge.iter())
                .map(|(&n, &m)| {
                    let total_len = n << protocol.num_vars;
                    let commitment_len = MOS::compute_comm_len(vk, total_len);
                    (
                        transcript.read_n_ec_points(commitment_len).unwrap(),
                        transcript.squeeze_n_challenges(m),
                    )
                })
                .unzip();

            (
                witnesses.into_iter().flatten().collect_vec(),
                challenges.into_iter().flatten().collect_vec(),
            )
        };

        let beta = transcript.squeeze_challenge();

        let lookup_count = transcript.read_ec_point().unwrap();

        let gamma = transcript.squeeze_challenge();

        let lookup_perm_intermediate = transcript.read_ec_point().unwrap();

        let alpha = transcript.squeeze_challenge();
        let y = transcript.squeeze_challenge();

        let constraint_degree = protocol.constraint_expression.degree();
        let constraint_sum_check_proof =
            SC::read_proof(protocol.num_vars, constraint_degree, transcript)?;

        let intermediate_evals = transcript.read_n_scalars(protocol.num_all_polys)?;
        let eta = transcript.squeeze_challenge();

        let opening_degree = protocol.opening_expression.degree();
        let opening_sum_check_proof =
            SC::read_proof(protocol.num_vars, opening_degree, transcript)?;

        let final_evals = transcript.read_n_scalars(protocol.num_all_polys)?;

        let eval_groups = reorder_into_groups(evals, &protocol.segment_groups);
        let pcs = MOS::read_proof(vk, protocol, &eval_groups, transcript)?;

        Ok(Self::Proof {
            witness_comms: witnesses,
            challenges,
            beta,
            lookup_count_comm: lookup_count,
            gamma,
            lookup_perm_intermediate_comm: lookup_perm_intermediate,
            alpha,
            y,
            constraint_sum_check_proof,
            intermediate_evals,
            eta,
            opening_sum_check_proof,
            final_evals,
            pcs,
        })
    }

    fn verify(
        vk: &Self::VerifyingKey,
        protocol: &Protocol<C, L>,
        instances: &[Vec<<L>::LoadedScalar>],
        proof: &Self::Proof,
    ) -> Result<Self::Output, Error> {
        let loader = proof.eta.loader();
        let constraint_degree = protocol.constraint_expression.degree();
        let x = SC::verify(
            &proof.constraint_sum_check_proof,
            &protocol.constraint_expression,
            &proof.intermediate_evals,
            &proof.challenges,
            &[&proof.y],
            loader.load_zero(),
            protocol.num_vars,
            constraint_degree,
        )?;
        let intermediate_evals = proof.intermediate_evals[protocol.num_instance.len()..];
        let powers_of_eta = proof.eta.powers(intermediate_evals.len());
        let random_combined_eval = loader.sum_products(
            &powers_of_eta
                .iter()
                .zip(intermediate_evals.iter().rev())
                .collect_vec(),
        );
        let opening_degree = protocol.opening_expression.degree();
        let z = SC::verify(
            &proof.opening_sum_check_proof,
            &protocol.opening_expression,
            &proof.final_evals,
            &proof.challenges,
            &[x],
            &random_combined_eval,
            &protocol.num_vars,
            opening_degree,
        )?;
        let comms = iter::empty()
            .chain(iter::once(&vk.prep_perm_comm))
            .chain(&proof.witness_comms)
            .chain(iter::once(&proof.lookup_count_comm))
            .chain(iter::once(&proof.lookup_perm_intermediate_comm))
            .collect_vec();
        MOS::verify(vk, &comms, &z, &proof.pcs)?;
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct HyperPlonkVerifyingKey<
    C: CurveAffine,
    L: Loader<C>,
    OS: OpenScheme<C, L>,
    MOS: MultiOpenScheme<C, L, OS>,
> {
    pub prep_perm_comm: L::LoadedEcPoint,
    pub pcs_vk: MOS::VerifyingKey,
}

#[derive(Clone, Debug)]
pub struct HyperPlonkProof<C, L, SRC, SC, OS, MOS>
where
    C: CurveAffine,
    L: Loader<C>,
    SRC: SumcheckRoundVerifier<C, L>,
    SC: SumcheckVerifier<C, L, SRC>,
    OS: OpenScheme<C, L>,
    MOS: MultiOpenScheme<C, L, OS>,
{
    pub witness_comms: Vec<L::LoadedEcPoint>,
    pub challenges: Vec<L::LoadedScalar>,
    pub beta: L::LoadedScalar,
    pub lookup_count_comm: L::LoadedScalar,
    pub gamma: L::LoadedScalar,
    pub lookup_perm_intermediate_comm: L::LoadedEcPoint,
    pub alpha: L::LoadedScalar,
    pub y: Vec<L::LoadedScalar>,
    pub constraint_sum_check_proof: SC::Proof,
    pub intermediate_evals: Vec<L::LoadedScalar>,
    pub eta: L::LoadedScalar,
    pub opening_sum_check_proof: SC::Proof,
    pub final_evals: Vec<L::LoadedScalar>,
    pub pcs: MOS::Proof,
}

fn langranges<C, L>(
    protocol: &Protocol<C, L>,
    instances: &[Vec<L::LoadedScalar>],
) -> impl IntoIterator<Item = i32>
where
    C: CurveAffine,
    L: Loader<C>,
{
    let instance_eval_lagrange = protocol.instance_committing_key.is_none().then(|| {
        let queries = {
            let offset = protocol.preprocessed.len();
            let range = offset..offset + protocol.num_instance.len();
            protocol
                .quotient
                .numerator
                .used_query()
                .into_iter()
                .filter(move |query| range.contains(&query.poly))
        };
        let (min_rotation, max_rotation) = queries.fold((0, 0), |(min, max), query| {
            if query.rotation.0 < min {
                (query.rotation.0, max)
            } else if query.rotation.0 > max {
                (min, query.rotation.0)
            } else {
                (min, max)
            }
        });
        let max_instance_len =
            Iterator::max(instances.iter().map(|instance| instance.len())).unwrap_or_default();
        -max_rotation..max_instance_len as i32 + min_rotation.abs()
    });
    protocol
        .quotient
        .numerator
        .used_langrange()
        .into_iter()
        .chain(instance_eval_lagrange.into_iter().flatten())
}

pub(super) fn reorder_into_groups<T: Clone>(
    arr: Vec<T>,
    segment_groups: &[Vec<(usize, usize)>], // offset and length
) -> Vec<Vec<T>> {
    let mut groups = vec![];
    for segments in segment_groups {
        let mut group = vec![];
        for (offset, length) in segments {
            group.extend(arr[*offset..*offset + *length].to_vec());
        }
        groups.push(group);
    }
    groups
}
