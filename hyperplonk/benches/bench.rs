use std::{fs::File, io, time::Instant};

use ark_bls12_381::{G1Affine as C, Fr};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Write};
use ark_std::test_rng;
use hyperplonk::{
    prelude::{CustomizedGates, HyperPlonkErrors, MockCircuit},
    HyperPlonkSNARK,
};
use subroutines::{pcs::PolynomialCommitmentScheme, poly_iop::PolyIOP, MultilinearHyraxPCS, HyraxSRS};

const SUPPORTED_SIZE: usize = 20;
const MIN_NUM_VARS: usize = 8;
const MAX_NUM_VARS: usize = 20;
const MIN_CUSTOM_DEGREE: usize = 1;
const MAX_CUSTOM_DEGREE: usize = 32;
const HIGH_DEGREE_TEST_NV: usize = 15;

fn main() -> Result<(), HyperPlonkErrors> {
    let thread = rayon::current_num_threads();
    println!("start benchmark with #{} threads", thread);
    let mut rng = test_rng();
    let pcs_srs = {
        match read_srs() {
            Ok(p) => p,
            Err(_e) => {
                let srs =
                    MultilinearHyraxPCS::<C>::gen_srs_for_testing(&mut rng, SUPPORTED_SIZE)?;
                write_srs(&srs);
                srs
            }
        }
    };
    bench_jellyfish_plonk(&pcs_srs, thread)?;
    println!();
    bench_vanilla_plonk(&pcs_srs, thread)?;
    println!();
    for degree in MIN_CUSTOM_DEGREE..=MAX_CUSTOM_DEGREE {
        bench_high_degree_plonk(&pcs_srs, degree, thread)?;
        println!();
    }
    println!();

    Ok(())
}

fn read_srs() -> Result<HyraxSRS<C>, io::Error> {
    let mut f = File::open("srs.params")?;
    Ok(HyraxSRS::<C>::deserialize_unchecked(&mut f).unwrap())
}

fn write_srs(pcs_srs: &HyraxSRS<C>) {
    let mut f = File::create("srs.params").unwrap();
    pcs_srs.serialize_uncompressed(&mut f).unwrap();
}

fn bench_vanilla_plonk(
    pcs_srs: &HyraxSRS<C>,
    thread: usize,
) -> Result<(), HyperPlonkErrors> {
    let filename = format!("vanilla threads {}.txt", thread);
    let mut file = File::create(filename).unwrap();
    for nv in MIN_NUM_VARS..=MAX_NUM_VARS {
        let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
        bench_mock_circuit_zkp_helper(&mut file, nv, &vanilla_gate, &pcs_srs)?;
    }

    Ok(())
}

fn bench_jellyfish_plonk(
    pcs_srs: &HyraxSRS<C>,
    thread: usize,
) -> Result<(), HyperPlonkErrors> {
    let filename = format!("jellyfish threads {}.txt", thread);
    let mut file = File::create(filename).unwrap();
    for nv in MIN_NUM_VARS..=MAX_NUM_VARS {
        let jf_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
        bench_mock_circuit_zkp_helper(&mut file, nv, &jf_gate, &pcs_srs)?;
    }

    Ok(())
}

fn bench_high_degree_plonk(
    pcs_srs: &HyraxSRS<C>,
    degree: usize,
    thread: usize,
) -> Result<(), HyperPlonkErrors> {
    let filename = format!("high degree {} thread {}.txt", degree, thread);
    let mut file = File::create(filename).unwrap();
    println!("custom gate of degree {}", degree);
    let vanilla_gate = CustomizedGates::mock_gate(2, degree);
    bench_mock_circuit_zkp_helper(&mut file, HIGH_DEGREE_TEST_NV, &vanilla_gate, &pcs_srs)?;

    Ok(())
}

fn bench_mock_circuit_zkp_helper(
    file: &mut File,
    nv: usize,
    gate: &CustomizedGates,
    pcs_srs: &HyraxSRS<C>,
) -> Result<(), HyperPlonkErrors> {
    let repetition = if nv < 10 {
        5
    } else if nv < 20 {
        2
    } else {
        1
    };

    //==========================================================
    let circuit = MockCircuit::<Fr>::new(1 << nv, gate);
    assert!(circuit.is_satisfied());
    let index = circuit.index;
    //==========================================================
    // generate pk and vks
    let start = Instant::now();
    for _ in 0..repetition {
        let (_pk, _vk) = <PolyIOP<Fr> as HyperPlonkSNARK<
            C,
            MultilinearHyraxPCS<C>,
        >>::preprocess(&index, &pcs_srs)?;
    }
    println!(
        "key extraction for {} variables: {} us",
        nv,
        start.elapsed().as_micros() / repetition as u128
    );
    let (pk, vk) =
        <PolyIOP<Fr> as HyperPlonkSNARK<C, MultilinearHyraxPCS<C>>>::preprocess(
            &index, &pcs_srs,
        )?;
    //==========================================================
    // generate a proof
    let start = Instant::now();
    for _ in 0..repetition {
        let _proof =
            <PolyIOP<Fr> as HyperPlonkSNARK<C, MultilinearHyraxPCS<C>>>::prove(
                &pk,
                &circuit.public_inputs,
                &circuit.witnesses,
            )?;
    }
    let t = start.elapsed().as_micros() / repetition as u128;
    println!(
        "proving for {} variables: {} us",
        nv,
        start.elapsed().as_micros() / repetition as u128
    );
    file.write_all(format!("{} {}\n", nv, t).as_ref()).unwrap();

    let proof = <PolyIOP<Fr> as HyperPlonkSNARK<C, MultilinearHyraxPCS<C>>>::prove(
        &pk,
        &circuit.public_inputs,
        &circuit.witnesses,
    )?;
    //==========================================================
    // verify a proof
    let start = Instant::now();
    for _ in 0..repetition {
        let verify =
            <PolyIOP<Fr> as HyperPlonkSNARK<C, MultilinearHyraxPCS<C>>>::verify(
                &vk,
                &circuit.public_inputs,
                &proof,
            )?;
        assert!(verify);
    }
    println!(
        "verifying for {} variables: {} us",
        nv,
        start.elapsed().as_micros() / repetition as u128
    );
    Ok(())
}
