use ark_ff::PrimeField;
use rand_chacha::rand_core::OsRng;
use transcript::IOPTranscript;

pub struct RandomTape<F: PrimeField> {
    tape: IOPTranscript<F>,
}

impl<F: PrimeField> RandomTape<F> {
    pub fn new(name: &'static [u8]) -> Self {
        let tape = {
            let mut csprng: OsRng = OsRng;
            let mut tape = IOPTranscript::new(name);
            tape.append_field_element(b"init_randomness", &F::rand(&mut csprng))
                .unwrap();
            tape
        };
        Self { tape }
    }

    pub fn random_scalar(&mut self, label: &'static [u8]) -> F {
        self.tape.get_and_append_challenge(label).unwrap()
    }

    pub fn random_vector(&mut self, label: &'static [u8], len: usize) -> Vec<F> {
        self.tape
            .get_and_append_challenge_vectors(label, len)
            .unwrap()
    }
}
