use std::{io, iter};

use halo2_proofs::{
    circuit,
    curves::CurveAffine,
    ff::{FromUniformBytes, PrimeField},
    plonk::{Any, ConstraintSystem},
    transcript::EncodedChallenge,
};
use hyperplonk::backend::{
    pcs::PolynomialCommitmentScheme,
    piop::HyperPlonkVerifierParam,
    util::{
        expression::{Query, Rotation},
        transcript::{FieldTranscript, Transcript},
    },
    PlonkishCircuitInfo,
};
use itertools::max;
use num_integer::div_ceil;

use crate::{
    halo2_verifier::{loader::Loader, util::protocol::Expression},
    hyperplonk_verifier::Protocol,
};

#[derive(Clone, Debug, Default)]
pub struct Config {
    pub zk: bool,
    pub query_instance: bool,
    pub num_proof: usize,
    pub num_instance: Vec<usize>,
}

impl Config {
    pub fn hyrax() -> Self {
        Self {
            zk: false, // because the implemented hyperplonk is not zero-knowledge
            query_instance: false,
            num_proof: 1,
            ..Default::default()
        }
    }

    pub fn set_zk(mut self, zk: bool) -> Self {
        self.zk = zk;
        self
    }

    pub fn set_query_instance(mut self, query_instance: bool) -> Self {
        self.query_instance = query_instance;
        self
    }

    pub fn with_num_proof(mut self, num_proof: usize) -> Self {
        assert!(num_proof > 0);
        self.num_proof = num_proof;
        self
    }

    pub fn with_num_instance(mut self, num_instance: Vec<usize>) -> Self {
        self.num_instance = num_instance;
        self
    }
}

pub fn compile<'a, C: CurveAffine, L: Loader<C>, Pcs>(
    circuit_info: &PlonkishCircuitInfo<C::Scalar>,
    params: &HyperPlonkVerifierParam<C::Scalar, Pcs>,
    config: Config,
) -> (Protocol<C, L>, HyperPlonkVerifyingKey<C, L>)
where
    Pcs: PolynomialCommitmentScheme<C::Scalar>,
{
    let cs = ConstraintSystem::default();
    let Config {
        zk,
        query_instance,
        num_proof,
        num_instance,
    } = config;
    let lookup_degree = max(circuit_info.lookups.iter().map(|lookup| {
        let max_input_degree = max(lookup.iter().map(|(input, _)| input.degree())).unwrap_or(0);
        let max_table_degree = max(lookup.iter().map(|(_, table)| table.degree())).unwrap_or(0);
        max_input_degree + max_table_degree
    }))
    .unwrap_or(0);
    let max_degree = iter::empty()
        .chain(circuit_info.constraints.iter())
        .map(Expression::degree)
        .chain(Some(lookup_degree + 1))
        .chain(circuit_info.max_degree)
        .chain(Some(2))
        .unwrap();
    let polynomials =
        &Polynomials::new(&cs, zk, query_instance, num_instance, num_proof, max_degree);
    let prep_perm_comm = params
        .prep_perm_comm
        .iter()
        .cloned()
        .map(Into::into)
        .collect_vec();

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

    (
        Protocol {
            num_instance: params.num_instances,
            num_witness: params.num_witness_polys,
            num_challenge: params.num_challenges,
            num_vars: params.num_vars,
            num_preprocess: params.num_preprocess,
            num_all_polys: circuit_info.permutation_prod_offset()
                + circuit_info.num_permutation_chunks(),
            max_degree,
            segment_groups,
            constraint_expression: polynomials.convert(&params.constraint_expression),
            opening_expression: polynomials.convert(&params.open_expression),
        },
        HyperPlonkVerifyingKey { prep_perm_comm },
    )
}

struct Polynomials<'a, F: PrimeField> {
    cs: &'a ConstraintSystem<F>,
    zk: bool,
    query_instance: bool,
    num_proof: usize,
    num_fixed: usize,
    num_permutation_fixed: usize,
    num_instance: Vec<usize>,
    num_advice: Vec<usize>,
    num_challenge: Vec<usize>,
    advice_index: Vec<usize>,
    challenge_index: Vec<usize>,
    num_lookup_count: usize,
    num_permutation_frac: usize,
    num_permutation_prod: usize,
    num_lookup_h: usize,
}

impl<'a, F: PrimeField> Polynomials<'a, F> {
    fn new(
        cs: &'a ConstraintSystem<F>,
        zk: bool,
        query_instance: bool,
        num_instance: Vec<usize>,
        num_proof: usize,
        degree: usize,
    ) -> Self {
        let permutation_chunk_size = degree - 1;
        let num_permutation_chunks = div_ceil(
            &cs.permutation().get_columns().len(),
            &permutation_chunk_size,
        );

        let num_phase = *cs.advice_column_phase().iter().max().unwrap_or(&0) as usize + 1;
        let remapping = |phase: Vec<u8>| {
            let num = phase.iter().fold(vec![0; num_phase], |mut num, phase| {
                num[*phase as usize] += 1;
                num
            });
            let index = phase
                .iter()
                .scan(vec![0; num_phase], |state, phase| {
                    let index = state[*phase as usize];
                    state[*phase as usize] += 1;
                    Some(index)
                })
                .collect::<Vec<_>>();
            (num, index)
        };

        let (num_advice, advice_index) = remapping(cs.advice_column_phase());
        let (num_challenge, challenge_index) = remapping(cs.challenge_phase());
        assert_eq!(num_advice.iter().sum::<usize>(), cs.num_advice_columns());
        assert_eq!(num_challenge.iter().sum::<usize>(), cs.num_challenges());
        Self {
            cs,
            zk,
            query_instance,
            num_proof,
            num_fixed: cs.num_fixed_columns(),
            num_permutation_fixed: cs.permutation().get_columns().len(),
            num_instance,
            num_advice,
            num_challenge,
            advice_index,
            challenge_index,
            num_lookup_count: cs.lookups().len(),
            num_permutation_frac: num_permutation_chunks,
            num_permutation_prod: num_permutation_chunks,
            num_lookup_h: cs.lookups().len(),
        }
    }

    fn num_preprocessed(&self) -> usize {
        self.num_fixed + self.num_permutation_fixed
    }

    fn num_instance(&self) -> Vec<usize> {
        iter::repeat(self.num_instance.clone())
            .take(self.num_proof)
            .flatten()
            .collect()
    }

    fn num_witness(&self) -> Vec<usize> {
        iter::empty()
            .chain(
                self.num_advice
                    .clone()
                    .iter()
                    .map(|num| self.num_proof * num),
            )
            .chain([
                self.num_proof * self.num_lookup_permuted,
                self.num_proof * (self.num_permutation_z + self.num_lookup_z) + self.zk as usize,
            ])
            .collect()
    }

    fn num_challenge(&self) -> Vec<usize> {
        let mut num_challenge = self.num_challenge.clone();
        *num_challenge.last_mut().unwrap() += 1; // theta
        iter::empty()
            .chain(num_challenge)
            .chain([
                1, // beta
                1, // gamma
                1, // alpha
                1, // eta
            ])
            .collect()
    }

    fn instance_offset(&self) -> usize {
        self.num_preprocessed()
    }

    fn witness_offset(&self) -> usize {
        self.instance_offset() + self.num_instance().len()
    }

    fn cs_witness_offset(&self) -> usize {
        self.witness_offset()
            + self
                .num_witness()
                .iter()
                .take(self.num_advice.len())
                .sum::<usize>()
    }

    fn query<C: Into<Any> + Copy, R: Into<Rotation>>(
        &self,
        column_type: C,
        mut column_index: usize,
        rotation: R,
        t: usize,
    ) -> Query {
        let offset = match column_type.into() {
            Any::Fixed => 0,
            Any::Instance => self.instance_offset() + t * self.num_instance.len(),
            Any::Advice(advice) => {
                column_index = self.advice_index[column_index];
                let phase_offset = self.num_proof
                    * self.num_advice[..advice.phase() as usize]
                        .iter()
                        .sum::<usize>();
                self.witness_offset() + phase_offset + t * self.num_advice[advice.phase() as usize]
            },
        };
        Query::new(offset + column_index, rotation.into())
    }

    fn instance_queries(&'a self, t: usize) -> impl IntoIterator<Item = Query> + 'a {
        self.query_instance
            .then(|| {
                self.cs
                    .instance_queries()
                    .iter()
                    .map(move |(column, rotation)| {
                        self.query(*column.column_type(), column.index(), *rotation, t)
                    })
            })
            .into_iter()
            .flatten()
    }

    fn advice_queries(&'a self, t: usize) -> impl IntoIterator<Item = Query> + 'a {
        self.cs
            .advice_queries()
            .iter()
            .map(move |(column, rotation)| {
                self.query(*column.column_type(), column.index(), *rotation, t)
            })
    }

    fn fixed_queries(&'a self) -> impl IntoIterator<Item = Query> + 'a {
        self.cs
            .fixed_queries()
            .iter()
            .map(move |(column, rotation)| {
                self.query(*column.column_type(), column.index(), *rotation, 0)
            })
    }

    fn permutation_fixed_queries(&'a self) -> impl IntoIterator<Item = Query> + 'a {
        (0..self.num_permutation_fixed).map(|i| Query::new(self.num_fixed + i, 0))
    }

    fn convert(
        &self,
        expression: &hyperplonk::backend::util::expression::Expression<F>,
    ) -> Expression<F> {
        expression.evaluate(
            &|scalar| Expression::Constant(scalar),
            &|idx| Expression::CommonPolynomial(idx),
            &|query| todo!(),
            &|challenge| {
                let phase_offset = self.num_challenge[..challenge.phase() as usize]
                    .iter()
                    .sum::<usize>();
                Expression::Challenge(phase_offset + self.challenge_index[challenge.index()])
            },
            &|a| -a,
            &|a, b| a + b,
            &|a, b| a * b,
            &|a, scalar| a * scalar,
        )
    }

    fn system_challenge_offset(&self) -> usize {
        let num_challenge = self.num_challenge();
        num_challenge[..num_challenge.len() - 4].iter().sum()
    }

    fn beta(&self) -> Expression<F> {
        Expression::Challenge(self.system_challenge_offset())
    }

    fn gamma(&self) -> Expression<F> {
        Expression::Challenge(self.system_challenge_offset() + 1)
    }

    fn alpha(&self) -> Expression<F> {
        Expression::Challenge(self.system_challenge_offset() + 2)
    }

    fn eta(&self) -> Expression<F> {
        Expression::Challenge(self.system_challenge_offset() + 3)
    }
}

struct MockChallenge;

impl<C: CurveAffine> EncodedChallenge<C> for MockChallenge {
    type Input = ();

    fn new(_: &Self::Input) -> Self {
        unreachable!()
    }

    fn get_scalar(&self) -> C::Scalar {
        unreachable!()
    }
}

#[derive(Default)]
struct MockTranscript<F: PrimeField>(F);

impl<C: CurveAffine> FieldTranscript<C::Scalar> for MockTranscript<C::Scalar> {
    fn squeeze_challenge(&mut self) -> MockChallenge {
        unreachable!();
    }
    fn common_field_element(&mut self, scalar: C::Scalar) -> io::Result<()> {
        self.0 = scalar;
        Ok(())
    }

    fn squeeze_challenges(&mut self, n: usize) -> Vec<C> {
        unreachable!();
    }
}
