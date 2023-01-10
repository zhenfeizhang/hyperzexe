use ark_bls12_381::{Fr, G1Affine as G1};
use ark_ff::UniformRand;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{sync::Arc, test_rng};
use std::time::Instant;
use subroutines::pcs::prelude::{MultilinearHyraxPCS, PCSError, PolynomialCommitmentScheme};
use transcript::IOPTranscript;

fn main() -> Result<(), PCSError> {
    bench_pcs()
}

fn bench_pcs() -> Result<(), PCSError> {
    let mut rng = test_rng();

    // normal polynomials
    let uni_params = MultilinearHyraxPCS::<G1>::gen_srs_for_testing(&mut rng, 24)?;

    for nv in 4..25 {
        let repetition = if nv < 10 {
            10
        } else if nv < 20 {
            5
        } else {
            2
        };

        let poly = Arc::new(DenseMultilinearExtension::rand(nv, &mut rng));
        let (ck, vk) = (uni_params.clone(), uni_params.clone());

        let point: Vec<_> = (0..nv).map(|_| Fr::rand(&mut rng)).collect();

        // commit
        let com = {
            let start = Instant::now();
            for _ in 0..repetition {
                let _commit = MultilinearHyraxPCS::commit(&ck, &poly)?;
            }

            println!(
                "Hyrax commit for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );

            MultilinearHyraxPCS::commit(&ck, &poly)?
        };

        // open
        let (proof, value) = {
            let mut transcript = IOPTranscript::<Fr>::new(b"test");
            let start = Instant::now();
            for _ in 0..repetition {
                let _open = MultilinearHyraxPCS::open(&ck, &poly, &point, &mut transcript)?;
            }

            println!(
                "Hyrax open for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
            let mut transcript = IOPTranscript::<Fr>::new(b"test");
            MultilinearHyraxPCS::open(&ck, &poly, &point, &mut transcript)?
        };

        // verify
        {
            let mut transcript = IOPTranscript::<Fr>::new(b"test");
            let start = Instant::now();
            for _ in 0..repetition {
                assert!(MultilinearHyraxPCS::verify(
                    &vk,
                    &com,
                    &point,
                    &value,
                    &proof,
                    &mut transcript,
                )?);
            }
            println!(
                "Hyrax verify for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
        }

        println!("====================================");
    }

    Ok(())
}
