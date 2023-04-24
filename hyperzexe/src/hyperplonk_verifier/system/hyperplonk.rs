use crate::{halo2_verifier::util::protocol::Expression, hyperplonk_verifier::Protocol};

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
            zk: true,
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

pub fn compile<'a, C: CurveAffine, P: Params<'a, C>>(
    params: &P,
    vk: &VerifyingKey<C>,
    config: Config,
) -> Protocol<C>
where
    C::ScalarExt: FromUniformBytes<64>,
{
    Protocol {
        num_instance: todo!(),
        num_witness: todo!(),
        num_challenge: todo!(),
        num_vars: todo!(),
        num_preprocess: todo!(),
        constraint_expression: todo!(),
        opening_expression: todo!(),
        prep_perm_comm: todo!(),
        evaluations: todo!(),
        transcript_initial_state: todo!(),
    }
}

impl From<poly::Rotation> for Rotation {
    fn from(rotation: poly::Rotation) -> Rotation {
        Rotation(rotation.0)
    }
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
    permutation_chunk_size: usize,
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
    ) -> Self {
        unimplemented!();
    }

    fn num_preprocessed(&self) -> usize {
        unimplemented!();
    }

    fn num_instance(&self) -> Vec<usize> {
        unimplemented!();
    }

    fn num_witness(&self) -> Vec<usize> {
        unimplemented!();
    }

    fn num_challenge(&self) -> Vec<usize> {
        unimplemented!();
    }

    pub fn num_poly(&self) -> usize {
        unimplemented!();
    }

    pub fn instance_offset(&self) -> usize {
        unimplemented!();
    }

    pub fn preprocess_offset(&self) -> usize {
        unimplemented!();
    }

    pub fn witness_offset(&self) -> usize {
        unimplemented!();
    }

    pub fn permutation_offset(&self) -> usize {
        unimplemented!();
    }

    pub fn lookup_count_offset(&self) -> usize {
        unimplemented!();
    }

    pub fn lookup_h_offset(&self) -> usize {
        unimplemented!();
    }

    pub fn permutation_frac_offset(&self) -> usize {
        unimplemented!();
    }

    pub fn permutation_prod_offset(&self) -> usize {
        unimplemented!();
    }

    pub fn permutation_p1_p2_offset(&self) -> usize {
        unimplemented!();
    }

    fn num_permutation_chunks(&self) -> usize {
        unimplemented!();
    }

    fn permutation_chunk_size(&self) -> usize {
        unimplemented!();
    }

    fn constraint_expression(&self) -> Expression<F> {
        unimplemented!();
    }

    fn opening_expression(&self) -> Expression<F> {
        unimplemented!();
    }

    fn system_challenge_offset(&self) -> usize {
        let num_challenge = self.num_challenge();
        num_challenge[..num_challenge.len() - 3].iter().sum()
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

impl<C: CurveAffine> Transcript<C, MockChallenge> for MockTranscript<C::Scalar> {
    fn squeeze_challenge(&mut self) -> MockChallenge {
        unreachable!()
    }

    fn common_point(&mut self, _: C) -> io::Result<()> {
        unreachable!()
    }

    fn common_scalar(&mut self, scalar: C::Scalar) -> io::Result<()> {
        self.0 = scalar;
        Ok(())
    }
}

fn transcript_initial_state<C: CurveAffine>(vk: &VerifyingKey<C>) -> C::Scalar
where
    C::Scalar: FromUniformBytes<64>,
{
    let mut transcript = MockTranscript::default();
    vk.hash_into(&mut transcript).unwrap();
    transcript.0
}
