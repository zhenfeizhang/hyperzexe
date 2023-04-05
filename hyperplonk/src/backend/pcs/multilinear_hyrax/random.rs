use halo2_curves::group::ff::PrimeField;
use merlin::Transcript;
use rand::rngs::OsRng;

use crate::backend::util::arithmetic::fe_from_bytes_le;

pub struct RandomTape<F: PrimeField> {
    tape: Transcript,
    Phantom: std::marker::PhantomData<F>,
}

impl<F: PrimeField> RandomTape<F> {
    pub fn new(name: &'static [u8]) -> Self {
        let tape = {
            let mut csprng: OsRng = OsRng;
            let mut tape = Transcript::new(name);
            let rand_bytes = F::random(&mut csprng).to_repr();
            tape.append_message(b"init_randomness", rand_bytes.as_ref());
            tape
        };
        Self {
            tape,
            Phantom: std::marker::PhantomData,
        }
    }

    pub fn random_scalar(&mut self, label: &'static [u8]) -> F {
        let mut buf = [0u8; 64];
        self.tape.challenge_bytes(label, &mut buf);
        fe_from_bytes_le(&buf)
    }

    pub fn random_vector(&mut self, label: &'static [u8], len: usize) -> Vec<F> {
        (0..len)
            .map(|_i| self.random_scalar(label))
            .collect::<Vec<F>>()
    }
}
