use std::{marker::PhantomData, sync::Arc};

use halo2_curves::group::ff::PrimeField;
use itertools::Itertools;

use crate::{
    backend::{
        poly::multilinear::MultilinearPolynomial,
        util::{
            end_timer,
            expression::{Expression, Query, Rotation},
            start_timer,
        },
        PlonkishCircuitInfo,
    },
    Error,
};

use self::commit::{compute_frac_poly, compute_p1_p2_poly, compute_product_poly};

mod commit;

pub(super) struct PermutationCheck<F>(PhantomData<F>);

impl<F: PrimeField> PermutationCheck<F> {
    /// generate the expressions for the virtual polynomial
    ///    prod(x) - p1(x) * p2(x) + alpha * [frac(x) * g1(x) * ... * gk(x) -
    /// f1(x)
    /// * ... * fk(x)] where p1(x) = (1-x1) * frac(x2, ..., xn, 0) + x1 *
    ///   prod(x2, ..., xn, 0), p2(x) = (1-x1) * frac(x2, ..., xn, 1) + x1 *
    ///   prod(x2, ..., xn, 1)
    /// Returns proof.
    pub(super) fn build_constraints(
        circuit_info: &PlonkishCircuitInfo<F>,
        beta: &Expression<F>,
        gamma: &Expression<F>,
    ) -> Vec<Expression<F>> {
        let permutation_polys = circuit_info.permutation_polys();
        let chunk_size = circuit_info.permutation_chunk_size();
        let num_chunks = circuit_info.num_permutation_chunks();
        let permutation_offset = circuit_info.permutation_offset();
        let frac_offset = circuit_info.permutation_frac_offset();
        let prod_offset = circuit_info.permutation_prod_offset();
        let p1_p2_offset = circuit_info.permutation_p1_p2_offset();
        let l_minus_2 = Expression::<F>::lagrange(-2);
        let one = Expression::<F>::Constant(F::one());
        let polys = permutation_polys
            .iter()
            .map(|idx| Query::new(*idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .collect_vec();
        let ids = (0..polys.len())
            .map(|idx| Expression::identity(idx))
            .collect_vec();
        let permutations = (permutation_offset..)
            .map(|idx| Query::new(idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .take(permutation_polys.len())
            .collect_vec();
        let fracs = (frac_offset..)
            .map(|idx| Query::new(idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .take(num_chunks)
            .collect_vec();
        let prods = (prod_offset..)
            .map(|idx| Query::new(idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .take(num_chunks)
            .collect_vec();
        let p1_p2 = (p1_p2_offset..)
            .step_by(2)
            .map(|idx| {
                (
                    Query::new(idx, Rotation::cur()),
                    Query::new(idx + 1, Rotation::cur()),
                )
            })
            .map(|(p1, p2)| {
                (
                    Expression::<F>::Polynomial(p1),
                    Expression::<F>::Polynomial(p2),
                )
            })
            .take(num_chunks)
            .collect_vec();
        polys
            .chunks(chunk_size)
            .zip(ids.chunks(chunk_size))
            .zip(permutations.chunks(chunk_size))
            .zip(fracs.iter())
            .zip(prods.iter())
            .zip(p1_p2.iter())
            .map(|(((((polys, ids), permutations), frac), prod), (p1, p2))| {
                let q1 = prod - p1 * p2;

                let g = polys
                    .iter()
                    .zip(permutations)
                    .map(|(poly, permutation)| poly + beta * permutation + gamma)
                    .product::<Expression<_>>();
                let f = polys
                    .iter()
                    .zip(ids)
                    .map(|(poly, id)| poly + beta * id + gamma)
                    .product::<Expression<_>>();

                let q2 = frac * g - f;
                let q3 = &l_minus_2 * (prod - &one);
                vec![q1, q2, q3]
            })
            .flatten()
            .collect_vec()
    }

    /// Generate polynomials in the permutation argument. Some need committing.
    /// Output:
    /// - frac_polys, where frac_poly(x) = f1(x)...fm(x) / (g1(x)...gm(x)) for
    ///   each chunk.
    /// - prod_polys, where prod_poly(x) stores the binary tree of the product
    ///   of frac_poly entries for each chunk.
    /// - p1_p2_polys (not committed),
    pub(super) fn commit(
        max_degree: usize,
        permutation_polys: &[(usize, Arc<MultilinearPolynomial<F>>)],
        polys: &[Arc<MultilinearPolynomial<F>>],
        beta: &F,
        gamma: &F,
    ) -> Result<
        (
            Vec<Arc<MultilinearPolynomial<F>>>,
            Vec<Arc<MultilinearPolynomial<F>>>,
            Vec<(Arc<MultilinearPolynomial<F>>, Arc<MultilinearPolynomial<F>>)>,
        ),
        Error,
    > {
        let start = start_timer(|| "PermutationCheck::commit");
        if permutation_polys.is_empty() {
            return Err(Error::InvalidSnark(
                "No permutation polynomials".to_string(),
            ));
        }
        let frac_polys = compute_frac_poly(max_degree, permutation_polys, polys, beta, gamma)?;
        let prod_polys = compute_product_poly(frac_polys.as_slice())?;
        let p1_p2_polys = frac_polys
            .iter()
            .zip(prod_polys.iter())
            .map(|(frac, prod)| compute_p1_p2_poly(frac, prod))
            .collect_vec();
        end_timer(start);
        Ok((frac_polys, prod_polys, p1_p2_polys))
    }
}

#[cfg(test)]
mod tests {
    use halo2_curves::bn256::Fr as F;

    use super::*;

    #[test]
    fn test_build_constraints() {
        let circuit_info = {
            let mut circuit_info = PlonkishCircuitInfo::<F> {
                k: 1,
                num_witness_polys: vec![5],
                permutations: vec![
                    vec![(0, 0), (1, 0)],
                    vec![(1, 1), (0, 1)],
                    vec![(2, 0), (3, 1)],
                    vec![(3, 0), (2, 1)],
                ],
                max_degree: Some(5),
                ..Default::default()
            };
            circuit_info.initialize_permutation_info(circuit_info.permutation_polys().as_slice());
            circuit_info
        };
        let beta = Expression::Constant(F::one());
        let gamma = Expression::Constant(F::one() + F::one());
        let permutation_constraints =
            PermutationCheck::build_constraints(&circuit_info, &beta, &gamma);
        let polys = (0..)
            .map(|idx| Query::new(idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .take(5)
            .collect_vec();
        let perms = (circuit_info.permutation_offset()..)
            .take(4)
            .map(|idx| Query::new(idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .collect_vec();
        let fracs = (circuit_info.permutation_frac_offset()..)
            .take(1)
            .map(|idx| Query::new(idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .collect_vec();
        let prods = (circuit_info.permutation_prod_offset()..)
            .take(1)
            .map(|idx| Query::new(idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .collect_vec();
        let p1s = (circuit_info.permutation_p1_p2_offset()..)
            .step_by(2)
            .take(1)
            .map(|idx| Query::new(idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .collect_vec();
        let p2s = (circuit_info.permutation_p1_p2_offset() + 1..)
            .step_by(2)
            .take(1)
            .map(|idx| Query::new(idx, Rotation::cur()))
            .map(Expression::<F>::Polynomial)
            .collect_vec();
        let identities = (0..).map(Expression::<F>::identity).take(4).collect_vec();
        let l_minus_2 = &Expression::<F>::lagrange(-2);
        let one = &Expression::<F>::Constant(F::one());
        let fxs = (1..4)
            .map(|i| &polys[i] + &beta * &identities[i] + &gamma)
            .fold(&polys[0] + &beta * &identities[0] + &gamma, |acc, x| {
                acc * x
            });
        let gxs = (1..4)
            .map(|i| &polys[i] + &beta * &perms[i] + &gamma)
            .fold(&polys[0] + &beta * &perms[0] + &gamma, |acc, x| acc * x);
        let expected_constraints = vec![
            &prods[0] - &p1s[0] * &p2s[0],
            &fracs[0] * &gxs - &fxs,
            l_minus_2 * (&prods[0] - one),
        ];
        assert_eq!(permutation_constraints, expected_constraints);
    }
}
