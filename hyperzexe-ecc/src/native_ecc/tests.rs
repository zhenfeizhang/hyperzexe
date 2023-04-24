#![allow(unused_assignments, unused_imports, unused_variables)]
use super::*;
use crate::halo2_proofs::{
    circuit::*,
    dev::MockProver,
    group::Group,
    halo2curves::bn256::{Fq, Fr, G1Affine, G2Affine, G1, G2},
    plonk::*,
};
use halo2_base::{
    gates::range::RangeStrategy,
    utils::{bigint_to_fe, value_to_option, PrimeField},
    ContextParams, SKIP_FIRST_PASS,
};
use num_bigint::{BigInt, RandBigInt};
use std::{marker::PhantomData, ops::Neg};

#[derive(Default)]
pub struct MyCircuit<F> {
    pub P: Option<G1Affine>,
    pub Q: Option<G1Affine>,
    pub _marker: PhantomData<F>,
}

const NUM_ADVICE: usize = 2;
const NUM_FIXED: usize = 2;

impl<F: PrimeField> Circuit<F> for MyCircuit<F> {
    type Config = FpConfig<F, Fq>;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self {
            P: None,
            Q: None,
            _marker: PhantomData,
        }
    }

    fn configure(meta: &mut ConstraintSystem<F>) -> Self::Config {
        unimplemented!();
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<F>,
    ) -> Result<(), Error> {
        unimplemented!();
    }
}

#[cfg(test)]
#[test]
fn test_ecc() {
    let k = 23;
    let mut rng = rand::thread_rng();

    let P = Some(G1Affine::random(&mut rng));
    let Q = Some(G1Affine::random(&mut rng));

    let circuit = MyCircuit::<Fr> {
        P,
        Q,
        _marker: PhantomData,
    };

    let prover = MockProver::run(k, &circuit, vec![]).unwrap();
    prover.assert_satisfied();
}

#[cfg(feature = "dev-graph")]
#[cfg(test)]
#[test]
fn plot_ecc() {
    let k = 10;
    use plotters::prelude::*;

    let root = BitMapBackend::new("layout.png", (512, 16384)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let root = root.titled("Ecc Layout", ("sans-serif", 60)).unwrap();

    let circuit = MyCircuit::<Fr>::default();

    halo2_proofs::dev::CircuitLayout::default()
        .render(k, &circuit, &root)
        .unwrap();
}
