use arithmetic::{identity_permutation_mles, VPAuxInfo, VirtualPolynomial};
use ark_bls12_381::{G1Affine as G1, Fr};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::test_rng;
use rand_chacha::rand_core::RngCore;
use std::{marker::PhantomData, sync::Arc, time::Instant};
use subroutines::{
    pcs::{prelude::MultilinearHyraxPCS, PolynomialCommitmentScheme},
    poly_iop::prelude::{
        PermutationCheck, PolyIOP, PolyIOPErrors, ProductCheck, SumCheck, ZeroCheck, LookupCheck,
    },
};

type Hyrax = MultilinearHyraxPCS<G1>;

fn main() -> Result<(), PolyIOPErrors> {
    bench_lookup_check()
    // println!("\n\n");
    // bench_permutation_check()?;
    // println!("\n\n");
    // bench_sum_check()?;
    // println!("\n\n");
    // bench_prod_check()?;
    // println!("\n\n");
    // bench_zero_check()
}


fn bench_lookup_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();

    for nv in 4..20 {
        let srs = Hyrax::gen_srs_for_testing(&mut rng, nv + 1)?;
        let (pcs_param, _) = Hyrax::trim(&srs, None, Some(nv + 1))?;

        let repetition = if nv < 10 {
            100
        } else if nv < 20 {
            50
        } else {
            10
        };

        let ws = Arc::new(DenseMultilinearExtension::rand(nv, &mut rng));
        let fsx = vec![Arc::new(DenseMultilinearExtension::from_evaluations_vec(nv,
            (0..(1<<nv)).map(|_| ws.evaluations[(rng.next_u32() % (1 << nv)) as usize]).collect::<Vec<_>>(),
        ))];

        let proof =
            {
                let start = Instant::now();
                let mut transcript =
                    <PolyIOP<Fr> as LookupCheck<G1, Hyrax>>::init_transcript();
                transcript.append_message(b"testing", b"initializing transcript for testing")?;

                let (proof, _q_x, _frac_poly) = <PolyIOP<Fr> as LookupCheck<
                    G1,
                    Hyrax,
                >>::prove(
                    &pcs_param, &fsx, &ws, &mut transcript
                )?;

                println!(
                    "Lookup check proving time for {} variables: {} ns",
                    nv,
                    start.elapsed().as_nanos() / repetition as u128
                );
                proof
            };

        {
            let poly_info = VPAuxInfo {
                max_degree: 3,
                num_variables: nv,
                phantom: PhantomData::default(),
            };

            let start = Instant::now();
            let mut transcript =
                <PolyIOP<Fr> as LookupCheck<G1, Hyrax>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let _perm_check_sum_claim = <PolyIOP<Fr> as LookupCheck<G1, Hyrax>>::verify(
                &proof,
                &poly_info,
                &mut transcript,
            )?;
            println!(
                "Lookup check verification time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
        }

        println!("====================================");
    }

    Ok(())
}

fn bench_sum_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();
    for degree in 2..4 {
        for nv in 4..25 {
            let repetition = if nv < 10 {
                100
            } else if nv < 20 {
                50
            } else {
                10
            };

            let (poly, asserted_sum) =
                VirtualPolynomial::rand(nv, (degree, degree + 1), 2, &mut rng)?;
            let poly_info = poly.aux_info.clone();
            let proof = {
                let start = Instant::now();
                for _ in 0..repetition {
                    let mut transcript = <PolyIOP<Fr> as SumCheck<Fr>>::init_transcript();
                    let _proof = <PolyIOP<Fr> as SumCheck<Fr>>::prove(&poly, &mut transcript)?;
                }

                println!(
                    "sum check proving time for {} variables and {} degree: {} ns",
                    nv,
                    degree,
                    start.elapsed().as_nanos() / repetition as u128
                );
                let mut transcript = <PolyIOP<Fr> as SumCheck<Fr>>::init_transcript();
                <PolyIOP<Fr> as SumCheck<Fr>>::prove(&poly, &mut transcript)?
            };

            {
                let start = Instant::now();

                for _ in 0..repetition {
                    let mut transcript = <PolyIOP<Fr> as SumCheck<Fr>>::init_transcript();
                    let _subclaim = <PolyIOP<Fr> as SumCheck<Fr>>::verify(
                        asserted_sum,
                        &proof,
                        &poly_info,
                        &mut transcript,
                    )?;
                }
                println!(
                    "sum check verification time for {} variables and {} degree: {} ns",
                    nv,
                    degree,
                    start.elapsed().as_nanos() / repetition as u128
                );
            }

            println!("====================================");
        }
    }
    Ok(())
}

fn bench_zero_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();
    for degree in 2..4 {
        for nv in 4..20 {
            let repetition = if nv < 10 {
                100
            } else if nv < 20 {
                50
            } else {
                10
            };

            let poly = VirtualPolynomial::rand_zero(nv, (degree, degree + 1), 2, &mut rng)?;
            let poly_info = poly.aux_info.clone();
            let proof = {
                let start = Instant::now();
                let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
                transcript.append_message(b"testing", b"initializing transcript for testing")?;
                let proof = <PolyIOP<Fr> as ZeroCheck<Fr>>::prove(&poly, &mut transcript)?;

                println!(
                    "zero check proving time for {} variables and {} degree: {} ns",
                    nv,
                    degree,
                    start.elapsed().as_nanos() / repetition as u128
                );
                proof
            };

            {
                let start = Instant::now();
                let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
                transcript.append_message(b"testing", b"initializing transcript for testing")?;
                let _zero_subclaim =
                    <PolyIOP<Fr> as ZeroCheck<Fr>>::verify(&proof, &poly_info, &mut transcript)?;
                println!(
                    "zero check verification time for {} variables and {} degree: {} ns",
                    nv,
                    degree,
                    start.elapsed().as_nanos() / repetition as u128
                );
            }

            println!("====================================");
        }
    }
    Ok(())
}

fn bench_permutation_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();

    for nv in 4..20 {
        let srs = Hyrax::gen_srs_for_testing(&mut rng, nv + 1)?;
        let (pcs_param, _) = Hyrax::trim(&srs, None, Some(nv + 1))?;

        let repetition = if nv < 10 {
            100
        } else if nv < 20 {
            50
        } else {
            10
        };

        let ws = vec![Arc::new(DenseMultilinearExtension::rand(nv, &mut rng))];

        // identity map
        let perms = identity_permutation_mles(nv, 1);

        let proof =
            {
                let start = Instant::now();
                let mut transcript =
                    <PolyIOP<Fr> as PermutationCheck<G1, Hyrax>>::init_transcript();
                transcript.append_message(b"testing", b"initializing transcript for testing")?;

                let (proof, _q_x, _frac_poly) = <PolyIOP<Fr> as PermutationCheck<
                    G1,
                    Hyrax,
                >>::prove(
                    &pcs_param, &ws, &ws, &perms, &mut transcript
                )?;

                println!(
                    "permutation check proving time for {} variables: {} ns",
                    nv,
                    start.elapsed().as_nanos() / repetition as u128
                );
                proof
            };

        {
            let poly_info = VPAuxInfo {
                max_degree: 2,
                num_variables: nv,
                phantom: PhantomData::default(),
            };

            let start = Instant::now();
            let mut transcript =
                <PolyIOP<Fr> as PermutationCheck<G1, Hyrax>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let _perm_check_sum_claim = <PolyIOP<Fr> as PermutationCheck<G1, Hyrax>>::verify(
                &proof,
                &poly_info,
                &mut transcript,
            )?;
            println!(
                "permutation check verification time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
        }

        println!("====================================");
    }

    Ok(())
}

fn bench_prod_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();

    for nv in 4..20 {
        let srs = Hyrax::gen_srs_for_testing(&mut rng, nv + 1)?;
        let (pcs_param, _) = Hyrax::trim(&srs, None, Some(nv + 1))?;

        let repetition = if nv < 10 {
            100
        } else if nv < 20 {
            50
        } else {
            10
        };

        let f: DenseMultilinearExtension<Fr> = DenseMultilinearExtension::rand(nv, &mut rng);
        let mut g = f.clone();
        g.evaluations.reverse();
        let fs = vec![Arc::new(f)];
        let gs = vec![Arc::new(g)];

        let proof = {
            let start = Instant::now();
            let mut transcript = <PolyIOP<Fr> as ProductCheck<G1, Hyrax>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;

            let (proof, _prod_x, _frac_poly) =
                <PolyIOP<Fr> as ProductCheck<G1, Hyrax>>::prove(
                    &pcs_param,
                    &fs,
                    &gs,
                    &mut transcript,
                )?;

            println!(
                "product check proving time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
            proof
        };

        {
            let poly_info = VPAuxInfo {
                max_degree: 2,
                num_variables: nv,
                phantom: PhantomData::default(),
            };

            let start = Instant::now();
            let mut transcript = <PolyIOP<Fr> as ProductCheck<G1, Hyrax>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let _perm_check_sum_claim = <PolyIOP<Fr> as ProductCheck<G1, Hyrax>>::verify(
                &proof,
                &poly_info,
                &mut transcript,
            )?;
            println!(
                "product check verification time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
        }

        println!("====================================");
    }

    Ok(())
}
