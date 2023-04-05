use crate::{
    backend::util::{expression::Expression, transcript::TranscriptWrite, Itertools},
    Error,
};
use halo2_curves::group::ff::Field;
use rand::RngCore;
use std::{collections::BTreeSet, fmt::Debug, iter};

use self::{pcs::PolynomialCommitmentScheme, util::arithmetic::div_ceil};

pub mod pcs;
pub mod piop;
pub mod poly;
pub mod util;

pub trait PlonkishBackend<F, Pcs>: Clone + Debug
where
    F: Field,
    Pcs: PolynomialCommitmentScheme<F>,
{
    type ProverParam: Debug;
    type VerifierParam: Debug;
    type Proof: Debug;

    fn setup(size: usize, rng: impl RngCore) -> Result<Pcs::SRS, Error>;

    fn preprocess(
        param: &Pcs::SRS,
        circuit_info: PlonkishCircuitInfo<F>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), Error>;

    fn prove(
        pp: &Self::ProverParam,
        instances: &[&[F]],
        circuit: &impl PlonkishCircuit<F>,
        transcript: &mut impl TranscriptWrite<Pcs::Curve, F>,
        rng: impl RngCore,
    ) -> Result<Self::Proof, Error>;

    fn verify(
        vp: &Self::VerifierParam,
        instances: &[&[F]],
        proof: &Self::Proof,
        transcript: &mut impl TranscriptWrite<Pcs::Curve, F>,
        rng: impl RngCore,
    ) -> Result<(), Error>;
}

#[derive(Clone, Debug, Default)]
pub struct PlonkishCircuitInfo<F> {
    /// 2^k is the size of the circuit
    pub k: usize,
    /// Number of instnace value in each instance polynomial.
    pub num_instances: Vec<usize>,
    /// Preprocessed polynomials, which has index starts with offset
    /// `num_instances.len()`.
    pub preprocess_polys: Vec<Vec<F>>,
    /// Number of witness polynoimal in each phase.
    /// Witness polynomial index starts with offset `num_instances.len()` +
    /// `preprocess_polys.len()`.
    pub num_witness_polys: Vec<usize>,
    /// Number of challenge in each phase.
    pub num_challenges: Vec<usize>,
    /// Constraints.
    pub constraints: Vec<Expression<F>>,
    /// Each item inside outer vector repesents an independent vector lookup,
    /// which contains vector of tuples representing the input and table
    /// respectively.
    pub lookups: Vec<Vec<(Expression<F>, Expression<F>)>>,
    /// Each item inside outer vector repesents an closed permutation cycle,
    /// which contains vetor of tuples representing the polynomial index and
    /// row respectively.
    pub permutations: Vec<Vec<(usize, usize)>>,
    /// Maximum degree of constraints
    pub max_degree: Option<usize>,
    /// Number of permutation polynomials, -1 means unknown.
    num_permutation_polys: usize,
    /// Number of permutation chunks, -1 means unknown.
    permutation_chunk_size: usize,
    /// Whether permutation info is initialized.
    pub permute_info_initialized: bool,
}

impl<F: Clone> PlonkishCircuitInfo<F> {
    pub fn is_well_formed(&self) -> bool {
        let num_challenges = self.num_challenges.iter().sum::<usize>();
        let polys = iter::empty()
            .chain(self.expressions().flat_map(Expression::used_poly))
            .chain(self.permutation_polys())
            .collect::<BTreeSet<_>>();
        let challenges = iter::empty()
            .chain(self.expressions().flat_map(Expression::used_challenge))
            .collect::<BTreeSet<_>>();
        let num_poly = self.num_poly();
        // Same amount of phases
        self.num_witness_polys.len() == self.num_challenges.len()
            // Every phase has some witness polys
            && !self.num_witness_polys.iter().any(|n| *n == 0)
            // Every phase except the last one has some challenges after the witness polys are committed
            && !self.num_challenges[..self.num_challenges.len() - 1].iter().any(|n| *n == 0)
            // Polynomial indices are in range
            && (polys.is_empty() || *polys.last().unwrap() < num_poly)
            // Challenge indices are in range
            && (challenges.is_empty() || *challenges.last().unwrap() < num_challenges)
            // Every constraint has degree less equal than `max_degree`
            && self
                .max_degree
                .map(|max_degree| {
                    !self
                        .constraints
                        .iter()
                        .any(|constraint| constraint.degree() > max_degree)
                })
                .unwrap_or(true)
    }

    pub fn initialize_permutation_info(&mut self, permutation_polys: &[usize]) {
        if self.permute_info_initialized {
            return;
        }
        self.num_permutation_polys = permutation_polys.len();

        let lookup_degree = self
            .lookups
            .iter()
            .map(|lookup| {
                let max_input_degree = lookup
                    .iter()
                    .map(|(input, _)| input.degree())
                    .max()
                    .unwrap_or(0);
                let max_table_degree = lookup
                    .iter()
                    .map(|(_, table)| table.degree())
                    .max()
                    .unwrap_or(0);
                max_input_degree + max_table_degree
            })
            .max()
            .unwrap_or(0);

        let max_degree = iter::empty()
            .chain(self.constraints.iter())
            .map(Expression::degree)
            .chain(Some(lookup_degree + 1))
            .chain(self.max_degree)
            .chain(Some(2))
            .max()
            .unwrap();
        self.permutation_chunk_size = max_degree - 1;
        self.permute_info_initialized = true;
        if self.num_permutation_chunks() > 1 {
            panic!("Permutation with chunks more than 1 hasn't been implemented yet");
        }
    }

    pub fn permutation_polys(&self) -> Vec<usize> {
        self.permutations
            .iter()
            .flat_map(|cycle| cycle.iter().map(|(poly, _)| *poly))
            .unique()
            .sorted()
            .collect()
    }

    pub fn expressions(&self) -> impl Iterator<Item = &Expression<F>> {
        iter::empty().chain(self.constraints.iter()).chain(
            self.lookups
                .iter()
                .flat_map(|lookup| lookup.iter().flat_map(|(input, table)| [input, table])),
        )
    }

    pub fn num_poly(&self) -> usize {
        self.num_instances.len()
            + self.preprocess_polys.len()
            + self.num_witness_polys.iter().sum::<usize>()
    }

    pub fn instance_offset(&self) -> usize {
        0
    }

    pub fn preprocess_offset(&self) -> usize {
        self.num_instances.len()
    }

    pub fn witness_offset(&self) -> usize {
        self.preprocess_offset() + self.preprocess_polys.len()
    }

    pub fn permutation_offset(&self) -> usize {
        self.num_poly()
    }

    pub fn lookup_count_offset(&self) -> usize {
        if !self.permute_info_initialized {
            panic!("Permutation info is not initialized");
        }
        self.permutation_offset() + self.num_permutation_polys
    }

    pub fn lookup_h_offset(&self) -> usize {
        self.lookup_count_offset() + self.lookups.len()
    }

    pub fn permutation_frac_offset(&self) -> usize {
        if !self.permute_info_initialized {
            panic!("Permutation info is not initialized");
        }
        self.lookup_h_offset() + self.lookups.len()
    }

    pub fn permutation_prod_offset(&self) -> usize {
        self.permutation_frac_offset() + self.num_permutation_chunks()
    }

    pub fn permutation_p1_p2_offset(&self) -> usize {
        self.permutation_prod_offset() + self.num_permutation_chunks()
    }

    fn num_permutation_chunks(&self) -> usize {
        div_ceil(self.num_permutation_polys, self.permutation_chunk_size)
    }

    fn permutation_chunk_size(&self) -> usize {
        self.permutation_chunk_size
    }
}

pub trait PlonkishCircuit<F> {
    fn synthesize(&self, round: usize, challenges: &[F]) -> Result<Vec<Vec<F>>, Error>;
}

#[cfg(any(test, feature = "benchmark"))]
mod test {
    use crate::{backend::PlonkishCircuit, Error};

    impl<F: Clone> PlonkishCircuit<F> for Vec<Vec<F>> {
        fn synthesize(&self, round: usize, challenges: &[F]) -> Result<Vec<Vec<F>>, Error> {
            assert!(round == 0 && challenges.is_empty());
            Ok(self.to_vec())
        }
    }
}
