#![allow(non_snake_case)]
use halo2_curves::CurveAffine;
use rand::rngs::OsRng;

use crate::{backend::util::arithmetic::variable_base_msm, Error};
use halo2_curves::group::{ff::Field, Curve, Group};

#[derive(Clone, Debug)]
pub struct MultiCommitGens<C: CurveAffine> {
    pub n: usize,
    pub G: Vec<C>,
    pub h: C,
}

impl<C: CurveAffine> MultiCommitGens<C> {
    pub fn new(n: usize, label: &[u8]) -> Result<Self, Error> {
        let gens = Self::sample_generators(n + 1, label);

        Ok(MultiCommitGens {
            n,
            G: gens[0..n].to_vec(),
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
                h: self.h.clone(),
            },
            MultiCommitGens {
                n: G2.len(),
                G: G2.to_vec(),
                h: self.h,
            },
        )
    }

    fn sample_generators(num_generators: usize, _label: &[u8]) -> Vec<C> {
        let mut rng = OsRng;
        let mut generators = Vec::with_capacity(num_generators);
        for _ in 0..num_generators {
            loop {
                let g = C::generator() * C::Scalar::random(&mut rng);
                if !bool::from(g.is_identity()) {
                    generators.push(g);
                    break;
                }
            }
        }

        let mut res = vec![C::identity(); num_generators];
        C::Curve::batch_normalize(&generators, &mut res);
        res
    }
}

pub fn commit_array<C: CurveAffine>(
    scalars: &[C::Scalar],
    blind: &C::Scalar,
    gens_n: &MultiCommitGens<C>,
) -> C::Curve {
    assert_eq!(gens_n.n, scalars.len());
    // change the type of self to &[<C::Scalar as PrimeField>::BigInt]
    variable_base_msm(scalars, &gens_n.G) + gens_n.h.mul(blind)
}

pub fn commit_element<C: CurveAffine>(
    scalar: &C::Scalar,
    blind: &C::Scalar,
    gens_1: &MultiCommitGens<C>,
) -> C::Curve {
    assert_eq!(gens_1.n, 1);
    // change the type of self to &[<C::Scalar as PrimeField>::BigInt]
    gens_1.G[0].mul(*scalar) + gens_1.h.mul(*blind)
}
