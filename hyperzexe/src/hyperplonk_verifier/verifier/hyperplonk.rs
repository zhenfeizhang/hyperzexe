use std::marker::PhantomData;

use halo2_proofs::curves::CurveAffine;

use crate::{hyperplonk_verifier::{Protocol, pcs::OpenScheme}, halo2_verifier::util::transcript::TranscriptRead};

use super::HyperPlonkVerifier;

pub struct HyperPlonk<PCS, SC>(PhantomData<(PCS, SC)>);

impl<C, L, PCS, SC> HyperPlonkVerifier<C, L, PCS, SC> for HyperPlonk<PCS, SC>
where
    C: CurveAffine,
    L: Loader<C>,
    PCS: OpenScheme<C, L>,
    SC:
{
    type Proof = HyperPlonkProof<C, L, PCS, SC>;

    fn read_proof<T>(
        vk: &PCS::VerifyingKey,
        protocol: &Protocol<C, L>,
        instances: &[Vec<L::LoadedScalar>],
        transcript: &mut T,
    ) -> Self::Proof
    where
        T: TranscriptRead<C, L>,
    {
        HyperPlonkProof::read::<T>(vk, protocol, instances, transcript)
    }

    fn verify(
        vk: &PCS::VerifyingKey,
        protocol: &Protocol<C, L>,
        instances: &[Vec<<L>::LoadedScalar>],
        proof: &Self::Proof,
    ) -> () {
        let loader = proof.eta.loader();
        let constraint_sum_check_msgs = proof.constraint_sum_check_msgs.iter().map(|x| x.as_slice()).collect::<Vec<_>>();
        SC::verify(
            protocol.num_vars,
            protocol.constraint_expression,
            instances,
            &constraint_sum_check_msgs,
            &proof.challenges,
            &proof.y,
            loader.load_zero(),
            &proof.intermediate_evals,
            &proof.x,
        );
        let intermediate_evals = proof.intermediate_evals[protocol.num_instance.len()..];
        let powers_of_eta = proof.eta.powers(intermediate_evals.len());
        let random_combined_eval = loader.sum_products(&powers_of_eta.iter().zip(intermediate_evals.iter().rev()).collect_vec());
        let opening_sum_check_msgs = proof.opening_sum_check_msgs.iter().map(|x| x.as_slice()).collect::<Vec<_>>();
        let final_evals = proof.pcs....
        SC::verify(
            protocol.num_vars,
            protocol.opening_expression,
            instances,
            &opening_sum_check_msgs,
            &proof.challenges,
            &proof.x,
            random_combined_eval,
            &final_evals,
            &proof.z,
        );
        PCS::verify(vk, &proof.pcs);
    }
}

#[derive(Clone, Debug)]
pub struct HyperPlonkProof<C, L, PCS>
where
    C: CurveAffine,
    L: Loader<C>,
{
    pub witnesses: Vec<L::LoadedEcPoint>,
    pub challenges: Vec<L::LoadedScalar>,
    pub lookup_count: L::LoadedScalar,
    pub lookup_perm_intermediate: L::LoadedEcPoint,
    pub beta: L::LoadedScalar,
    pub gamma: L::LoadedScalar,
    pub alpha: L::LoadedScalar,
    pub y: Vec<L::LoadedScalar>,
    pub constraint_sum_check_msgs: Vec<Vec<L::LoadedScalar>>,
    pub x: Vec<L::LoadedScalar>,
    pub intermediate_evals: Vec<L::LoadedScalar>,
    pub eta: L::LoadedScalar,
    pub opening_sum_check_msgs: Vec<Vec<L::LoadedScalar>>,
    pub z: L::LoadedScalar,
    pub pcs: PCS::BatchProof,
}

impl<C, L, PCS> HyperPlonkProof<C, L, PCS>
where
    PCS: OpenScheme<C, L>,
    C: CurveAffine,
    L: Loader<C>,
{
    pub fn read<T>(
        svk: &PCS::VerifyingKey,
        protocol: &Protocol<C, L>,
        instances: &[Vec<L::LoadedScalar>],
        transcript: &mut T,
    ) -> Self
    where
        T: TranscriptRead<C, L>,
    {
        if let Some(transcript_initial_state) = &protocol.transcript_initial_state {
            transcript.common_scalar(transcript_initial_state).unwrap();
        }

        debug_assert_eq!(
            protocol.num_instance,
            instances.iter().map(|instances| instances.len()).collect_vec(),
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
                    (transcript.read_n_ec_points(compute_comm_len(n)).unwrap(), transcript.squeeze_n_challenges(m))
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

        let (constraint_sum_check_msgs, x) = SC::read(transcript);

        // let constraint_sum_check_msgs = Vec::new();

        // let mut x = Vec::with_capacity(protocol.num_vars);
        // for i in 0..protocol.num_vars {
        //     constraint_sum_check_msgs.push(transcript.read_n_scalars(protocol.constraint_expression.degree())?);
        //     x.push(transcript.squeeze_challenge());
        // }

        let intermediate_evals = transcript.read_n_scalars(protocol.num_vars.num_intermediate_evals)?;
        let eta = transcript.squeeze_challenge();

        let (opening_sum_check_msgs, z) = SC::read(transcript);
        // let mut z = Vec::with_capacity(protocol.num_vars);
        // let mut opening_sum_check_msgs = Vec::with_capacity(protocol.num_vars);
        // for i in 0..protocol.num_vars {
        //     opening_sum_check_msgs.push(transcript.read_n_scalars(protocol.opening_expression.degree())?);
        //     z.push(transcript.squeeze_challenge());
        // }

        let evals = Vec::new();
        evals.extend(vec![loader.load_zero(); protocol.num_instance]);

        let segment_groups = {
            let mut offset = witness_offset;
            let mut witness_offsets = Vec::new();
            for num_witness_poly in pp.num_witness_polys.iter() {
                witness_offsets.push(vec![(offset, *num_witness_poly)]);
                offset += *num_witness_poly;
            }
            let mut res = vec![
                vec![
                    (preprocess_offset, pp.preprocess_polys.len()),
                    (permutation_offset, permutation_polys.len()),
                ]];
            res.extend(witness_offsets);
            res.extend(vec![
                vec![(lookup_count_offset, lookup_count_polys.len())],
                vec![
                    (lookup_h_offset, lookup_h_polys.len()),
                    (permutation_frac_offset, permutation_frac_polys.len()),
                    (permutation_prod_offset, permutation_prod_polys.len()),
                ],
            ]);
            res
        };

        let eval_groups = reorder_into_groups(evals, &segment_groups);
        let pcs = PCS::read(svk, protocol, transcript, &eval_groups)?;

        Ok(Self {
            witnesses,
            challenges,
            lookup_count,
            lookup_perm_intermediate,
            beta,
            gamma,
            alpha,
            y,
            constraint_sum_check_msgs,
            x,
            intermediate_evals,
            eta,
            opening_sum_check_msgs,
            z,
            pcs,
        })
    }

    pub fn empty_queries(protocol: &Protocol<C, L>) -> Vec<pcs::Query<C::Scalar>> {
        protocol
            .queries
            .iter()
            .map(|query| pcs::Query {
                poly: query.poly,
                shift: protocol.domain.rotate_scalar(C::Scalar::ONE, query.rotation),
                eval: (),
            })
            .collect()
    }

    fn queries(
        &self,
        protocol: &Protocol<C, L>,
        mut evaluations: FxHashMap<Query, L::LoadedScalar>,
    ) -> Vec<pcs::Query<C::Scalar, L::LoadedScalar>> {
        Self::empty_queries(protocol)
            .into_iter()
            .zip(protocol.queries.iter().map(|query| evaluations.remove(query).unwrap()))
            .map(|(query, eval)| query.with_evaluation(eval))
            .collect()
    }

    fn commitments<'a>(
        &'a self,
        protocol: &'a Protocol<C, L>,
        common_poly_eval: &CommonPolynomialEvaluation<C, L>,
        evaluations: &mut FxHashMap<Query, L::LoadedScalar>,
    ) -> Vec<Msm<C, L>> {
        let loader = common_poly_eval.zn().loader();
        let mut commitments = iter::empty()
            .chain(protocol.preprocessed.iter().map(Msm::base))
            .chain(
                self.committed_instances
                    .as_ref()
                    .map(|committed_instances| {
                        committed_instances.iter().map(Msm::base).collect_vec()
                    })
                    .unwrap_or_else(|| {
                        iter::repeat_with(Default::default)
                            .take(protocol.num_instance.len())
                            .collect_vec()
                    }),
            )
            .chain(self.witnesses.iter().map(Msm::base))
            .collect_vec();

        let numerator = protocol.quotient.numerator.evaluate(
            &|scalar| Msm::constant(loader.load_const(&scalar)),
            &|poly| Msm::constant(common_poly_eval.get(poly).clone()),
            &|query| {
                evaluations
                    .get(&query)
                    .cloned()
                    .map(Msm::constant)
                    .or_else(|| {
                        (query.rotation == Rotation::cur())
                            .then(|| commitments.get(query.poly).cloned())
                            .flatten()
                    })
                    .ok_or(Error::InvalidQuery(query))
                    .unwrap()
            },
            &|index| {
                self.challenges
                    .get(index)
                    .cloned()
                    .map(Msm::constant)
                    .ok_or(Error::InvalidChallenge(index))
                    .unwrap()
            },
            &|a| -a,
            &|a, b| a + b,
            &|a, b| match (a.size(), b.size()) {
                (0, _) => b * &a.try_into_constant().unwrap(),
                (_, 0) => a * &b.try_into_constant().unwrap(),
                (_, _) => panic!("{:?}", Error::InvalidLinearization),
            },
            &|a, scalar| a * &loader.load_const(&scalar),
        );

        let quotient_query = Query::new(
            protocol.preprocessed.len() + protocol.num_instance.len() + self.witnesses.len(),
            Rotation::cur(),
        );
        let quotient = common_poly_eval
            .zn()
            .pow_const(protocol.quotient.chunk_degree as u64)
            .powers(self.quotients.len())
            .into_iter()
            .zip(self.quotients.iter().map(Msm::base))
            .map(|(coeff, chunk)| chunk * &coeff)
            .sum::<Msm<_, _>>();
        match protocol.linearization {
            Some(LinearizationStrategy::WithoutConstant) => {
                let linearization_query = Query::new(quotient_query.poly + 1, Rotation::cur());
                let (msm, constant) = numerator.split();
                commitments.push(quotient);
                commitments.push(msm);
                evaluations.insert(
                    quotient_query,
                    (constant.unwrap_or_else(|| loader.load_zero())
                        + evaluations.get(&linearization_query).unwrap())
                        * common_poly_eval.zn_minus_one_inv(),
                );
            }
            Some(LinearizationStrategy::MinusVanishingTimesQuotient) => {
                let (msm, constant) =
                    (numerator - quotient * common_poly_eval.zn_minus_one()).split();
                commitments.push(msm);
                evaluations.insert(quotient_query, constant.unwrap_or_else(|| loader.load_zero()));
            }
            None => {
                commitments.push(quotient);
                evaluations.insert(
                    quotient_query,
                    numerator.try_into_constant().ok_or(Error::InvalidLinearization).unwrap()
                        * common_poly_eval.zn_minus_one_inv(),
                );
            }
        }

        commitments
    }

    fn evaluations(
        &self,
        protocol: &Protocol<C, L>,
        instances: &[Vec<L::LoadedScalar>],
        common_poly_eval: &CommonPolynomialEvaluation<C, L>,
    ) -> FxHashMap<Query, L::LoadedScalar> {
        let loader = common_poly_eval.zn().loader();
        let instance_evals = protocol.instance_committing_key.is_none().then(|| {
            let offset = protocol.preprocessed.len();
            let queries = {
                let range = offset..offset + protocol.num_instance.len();
                protocol
                    .quotient
                    .numerator
                    .used_query()
                    .into_iter()
                    .filter(move |query| range.contains(&query.poly))
            };
            queries
                .map(move |query| {
                    let instances = instances[query.poly - offset].iter();
                    let l_i_minus_r = (-query.rotation.0..)
                        .map(|i_minus_r| common_poly_eval.get(Lagrange(i_minus_r)));
                    let eval = loader.sum_products(&instances.zip(l_i_minus_r).collect_vec());
                    (query, eval)
                })
                .collect_vec()
        });

        iter::empty()
            .chain(instance_evals.into_iter().flatten())
            .chain(protocol.evaluations.iter().cloned().zip(self.evaluations.iter().cloned()))
            .collect()
    }
}

impl<C, MOS> CostEstimation<(C, MOS)> for Plonk<MOS>
where
    C: CurveAffine,
    MOS: MultiOpenScheme<C, NativeLoader> + CostEstimation<C, Input = Vec<pcs::Query<C::Scalar>>>,
{
    type Input = Protocol<C>;

    fn estimate_cost(protocol: &Protocol<C>) -> Cost {
        let plonk_cost = {
            let num_accumulator = protocol.accumulator_indices.len();
            let num_instance = protocol.num_instance.iter().sum();
            let num_commitment =
                protocol.num_witness.iter().sum::<usize>() + protocol.quotient.num_chunk();
            let num_evaluation = protocol.evaluations.len();
            let num_msm = protocol.preprocessed.len() + num_commitment + 1 + 2 * num_accumulator;
            Cost::new(num_instance, num_commitment, num_evaluation, num_msm)
        };
        let pcs_cost = {
            let queries = PlonkProof::<C, NativeLoader, MOS>::empty_queries(protocol);
            MOS::estimate_cost(&queries)
        };
        plonk_cost + pcs_cost
    }
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
