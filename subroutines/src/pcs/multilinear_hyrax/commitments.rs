#![allow(non_snake_case)]
use ark_ec::{msm::VariableBaseMSM, AffineCurve, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError, Write};

use blake2::Blake2s256;
use digest::Digest;
use std::io::Read;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::PCSError;

#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct MultiCommitGens<C: AffineCurve> {
    pub n: usize,
    pub G: Vec<C>,
    pub h: C,
}

impl<C: AffineCurve> MultiCommitGens<C> {
    pub fn new(n: usize, label: &[u8]) -> Result<Self, PCSError> {
        let gens = Self::sample_generators(n + 1, label);

        Ok(MultiCommitGens {
            n,
            G: gens[..n].to_vec(),
            h: gens[n],
        })
    }

    pub fn clone(&self) -> MultiCommitGens<C> {
        MultiCommitGens {
            n: self.n,
            h: self.h,
            G: self.G.clone(),
        }
    }

    pub fn split_at(&self, mid: usize) -> (MultiCommitGens<C>, MultiCommitGens<C>) {
        let (G1, G2) = self.G.split_at(mid);

        (
            MultiCommitGens {
                n: G1.len(),
                G: G1.to_vec(),
                h: self.h,
            },
            MultiCommitGens {
                n: G2.len(),
                G: G2.to_vec(),
                h: self.h,
            },
        )
    }

    // from https://github.com/arkworks-rs/poly-commit
    fn sample_generators(num_generators: usize, label: &[u8]) -> Vec<C> {
        let generators: Vec<_> = ark_std::cfg_into_iter!(0..num_generators)
            .map(|i| {
                let i = i as u64;
                let mut hash = Blake2s256::digest([label, &i.to_le_bytes()].concat().as_slice());
                let mut g = C::from_random_bytes(&hash);
                let mut j = 0u64;
                while g.is_none() {
                    // PROTOCOL NAME, i, j
                    let mut bytes = label.to_vec();
                    bytes.extend(i.to_le_bytes());
                    bytes.extend(j.to_le_bytes());
                    hash = Blake2s256::digest(bytes.as_slice());
                    g = C::from_random_bytes(&hash);
                    j += 1;
                }
                let generator = g.unwrap();
                generator.mul_by_cofactor_to_projective()
            })
            .collect();

        C::Projective::batch_normalization_into_affine(&generators)
    }
}

pub fn commit_array<C: AffineCurve>(
    scalars: &[C::ScalarField],
    blind: &C::ScalarField,
    gens_n: &MultiCommitGens<C>,
) -> C::Projective {
    assert_eq!(gens_n.n, scalars.len());
    // change the type of self to &[<C::ScalarField as PrimeField>::BigInt]
    let scalars: Vec<_> = scalars.into_iter().map(|s| s.into_repr()).collect();
    VariableBaseMSM::multi_scalar_mul(&gens_n.G, scalars.as_slice())
        + gens_n.h.mul(blind.into_repr())
}

pub fn commit_element<C: AffineCurve>(
    scalar: &C::ScalarField,
    blind: &C::ScalarField,
    gens_1: &MultiCommitGens<C>,
) -> C::Projective {
    assert_eq!(gens_1.n, 1);
    // change the type of self to &[<C::ScalarField as PrimeField>::BigInt]
    gens_1.G[0].mul(*scalar) + gens_1.h.mul(*blind)
}
