#![allow(non_snake_case)]
use crate::halo2_proofs::group::Curve;
use ark_std::{end_timer, start_timer};
use serde::{Deserialize, Serialize};
use std::{io::Write, marker::PhantomData};

use crate::{
    ecc::EccChip,
    fields::fp::FpStrategy,
    halo2_proofs::{
        circuit::{Layouter, SimpleFloorPlanner, Value},
        dev::MockProver,
        halo2curves::bn256::{pairing, Bn256, Fr, G1Affine},
        plonk::*,
        poly::{
            commitment::{Params, ParamsProver},
            kzg::{
                commitment::{KZGCommitmentScheme, ParamsKZG},
                multiopen::{ProverSHPLONK, VerifierSHPLONK},
                strategy::SingleStrategy,
            },
        },
        transcript::{
            Blake2bRead, Blake2bWrite, Challenge255, TranscriptReadBuffer, TranscriptWriterBuffer,
        },
    },
};
use halo2_base::{
    gates::GateInstructions,
    utils::{biguint_to_fe, fe_to_biguint, value_to_option, PrimeField},
    QuantumCell::Witness,
};
use num_bigint::BigUint;
use num_traits::Num;

pub mod ec_add;
pub mod fixed_base_msm;
pub mod msm;
