use core::hash::Hash;
use std::{collections::HashSet, sync::Arc};

use halo2_curves::group::ff::PrimeField;
use itertools::Itertools;

use crate::{
    backend::{
        poly::multilinear::MultilinearPolynomial,
        util::{
            arithmetic::{powers, BooleanHypercube},
            end_timer,
            expression::{CommonPolynomial, Expression, Query, Rotation},
            parallel::parallelize,
            start_timer,
        },
        PlonkishCircuitInfo,
    },
    Error,
};

use self::commit::{compute_count_in_input, compute_h};

mod commit;

#[derive(Clone, Debug, Default)]
pub(super) struct LookupCheck<F> {
    lookup_evals: Vec<(Arc<MultilinearPolynomial<F>>, Arc<MultilinearPolynomial<F>>)>,
}

impl<F: PrimeField + Hash> LookupCheck<F> {
    /// Generate the expression for the virtual polynomial
    ///    q1 = (h(x)(gx(x) + gamma) + count(x))\prod_{fx \in fxs}(fx(x) +
    /// gamma)
    ///       - gx(x)\sum_(fx \in fxs) \prod_{fx' \neq fx}(fx'(x) + gamma)
    /// Returns proof.
    /// - input: denominators = [f1(x) + gamma, .., fk(x) + gamma, g(x) +
    ///   gamma], h = \sum_{fx \in fxs} 1 / (fx + gamma) - count / (gx + beta)
    ///
    /// Cost: O(N)
    pub(super) fn build_constraints(
        circuit_info: &PlonkishCircuitInfo<F>,
        beta: &Expression<F>,
        gamma: &Expression<F>,
    ) -> (Vec<Expression<F>>, Vec<Expression<F>>) {
        let count_offset = circuit_info.lookup_count_offset();
        let h_offset = circuit_info.lookup_h_offset();
        let constraints = circuit_info
            .lookups
            .iter()
            .zip(count_offset..)
            .zip(h_offset..)
            .flat_map(|((lookup, count), h)| {
                let (inputs, tables) = lookup
                    .iter()
                    .map(|(input, table)| (input, table))
                    .unzip::<_, _, Vec<_>, Vec<_>>();
                let input = Expression::distribute_powers(inputs, beta);
                let table = Expression::distribute_powers(tables, beta);
                let f = input + gamma;
                let g = table + gamma;
                let h: Expression<F> = Expression::Polynomial(Query::new(h, Rotation::cur()));
                let count: Expression<F> =
                    Expression::Polynomial(Query::new(count, Rotation::cur()));
                // h = 1 / f - count / g
                vec![(h * &g + count) * f - &g]
            })
            .collect_vec();
        let sumchecks = (h_offset..(h_offset + circuit_info.lookups.len()))
            .map(|h| Expression::Polynomial(Query::new(h, Rotation::cur())))
            .collect_vec();
        (constraints, sumchecks)
    }

    pub(super) fn commit_counts(
        &mut self,
        lookups: &[&[(Expression<F>, Expression<F>)]],
        polys: &[Arc<MultilinearPolynomial<F>>],
        challenges: &[F],
        beta: &F,
    ) -> Result<Vec<Arc<MultilinearPolynomial<F>>>, Error> {
        if lookups.is_empty() {
            return Ok(vec![]);
        }
        let start = start_timer(|| "Lookup compute counts");
        if polys.is_empty() {
            return Err(Error::InvalidSnark("No polynomials".to_string()));
        }
        let num_vars = polys[0].num_vars();

        let bh = BooleanHypercube::new(num_vars);
        let expression = lookups
            .iter()
            .flat_map(|lookup| lookup.iter().map(|(input, table)| (input + table)))
            .sum::<Expression<_>>();
        let lagranges = {
            let bh = BooleanHypercube::new(num_vars).iter().collect_vec();
            expression
                .used_langrange()
                .into_iter()
                .map(|i| (i, bh[i.rem_euclid(1 << num_vars) as usize]))
                .collect::<HashSet<_>>()
        };
        let identities = {
            let max_used_identity = expression
                .used_identity()
                .into_iter()
                .max()
                .unwrap_or_default();
            (0..=max_used_identity)
                .map(|idx| (idx as u64) << num_vars)
                .collect_vec()
        };

        let compress = |powers_of_beta: &[F], expressions: &[&Expression<F>]| {
            let res = powers_of_beta
                .iter()
                .rev()
                .copied()
                .zip(expressions.iter().map(|expression| {
                    let mut compressed = vec![F::zero(); 1 << num_vars];
                    parallelize(&mut compressed, |(compressed, start)| {
                        for (b, compressed) in (start..).zip(compressed) {
                            *compressed = expression.evaluate(
                                &|constant| constant,
                                &|common_poly| match common_poly {
                                    CommonPolynomial::Lagrange(i) => lagranges
                                        .contains(&(i, b))
                                        .then(F::one)
                                        .unwrap_or_else(F::zero),
                                    CommonPolynomial::Identity(idx) => {
                                        F::from(b as u64 + identities[idx])
                                    },
                                    CommonPolynomial::EqXY(_) => unreachable!(),
                                    CommonPolynomial::EqXYShifted00(_) => unreachable!(),
                                    CommonPolynomial::EqXYShifted01(_) => unreachable!(),
                                    CommonPolynomial::EqXYShifted10(_) => unreachable!(),
                                    CommonPolynomial::EqXYShifted11(_) => unreachable!(),
                                },
                                &|query| polys[query.poly()][bh.rotate(b, query.rotation())],
                                &|challenge| challenges[challenge],
                                &|value| -value,
                                &|lhs, rhs| lhs + &rhs,
                                &|lhs, rhs| lhs * &rhs,
                                &|value, scalar| value * &scalar,
                            );
                        }
                    });
                    MultilinearPolynomial::new(compressed)
                }))
                .sum::<MultilinearPolynomial<_>>();
            Arc::new(res)
        };

        let counts = lookups
            .iter()
            .map(|lookup| {
                let powers_of_beta = powers(*beta).take(lookup.len()).collect_vec();
                let (inputs, tables) = lookup
                    .iter()
                    .map(|(input, table)| (input, table))
                    .unzip::<_, _, Vec<_>, Vec<_>>();
                let compressed_input_poly = compress(powers_of_beta.as_slice(), &inputs);
                let compressed_table_poly = compress(powers_of_beta.as_slice(), &tables);
                let count_poly = compute_count_in_input(
                    &[compressed_input_poly.clone()],
                    &compressed_table_poly,
                )?;
                self.lookup_evals
                    .push((compressed_input_poly, compressed_table_poly));
                Ok(count_poly)
            })
            .collect::<Result<Vec<_>, _>>()?;
        end_timer(start);
        Ok(counts)
    }

    pub(super) fn commit_hs(
        &mut self,
        count_polys: &[Arc<MultilinearPolynomial<F>>],
        gamma: &F,
    ) -> Result<Vec<Arc<MultilinearPolynomial<F>>>, Error> {
        let start = start_timer(|| "Lookup compute h polys");
        let mut h_polys = Vec::with_capacity(self.lookup_evals.len());
        for ((input_evals, table_evals), count_poly) in self.lookup_evals.iter().zip(count_polys) {
            let h_poly = compute_h(&[input_evals.clone()], table_evals, count_poly, gamma)?;
            h_polys.push(h_poly);
        }
        end_timer(start);
        Ok(h_polys)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::{
        poly::multilinear::arc_mpoly,
        util::{
            arithmetic::sum,
            expression::{Expression, Query, Rotation},
        },
    };
    use halo2_curves::{bn256::Fr as F, group::ff::Field};
    use rand::rngs::OsRng;
    #[test]
    fn test_commit_counts() {
        let tables = vec![
            arc_mpoly!(F, 1, 2, 3, 4),
            arc_mpoly!(F, 2, 4, 6, 8),
            arc_mpoly!(F, 3, 6, 9, 1),
        ];

        let inputs = vec![
            arc_mpoly!(F, 2, 3, 1, 1),
            arc_mpoly!(F, 4, 6, 2, 2),
            arc_mpoly!(F, 6, 9, 3, 3),
        ];

        let polys = [inputs, tables].concat();
        let poly_exps = (0..polys.len())
            .map(|poly| Expression::Polynomial(Query::new(poly, Rotation::cur())))
            .collect_vec();
        let two = Expression::Constant(F::one() + &F::one());
        let mut lookups = vec![
            vec![
                (&poly_exps[0] + &poly_exps[2], &poly_exps[3] + &poly_exps[5]),
                (
                    &poly_exps[1] * &two * &poly_exps[2],
                    &poly_exps[4] * &two * &poly_exps[5],
                ),
            ],
            vec![
                (&poly_exps[0] + &poly_exps[2], &poly_exps[3] + &poly_exps[5]),
                (
                    &poly_exps[1] * &two * &poly_exps[2],
                    &poly_exps[4] * &two * &poly_exps[5],
                ),
            ],
        ];
        // Right path
        let mut lookup_check = LookupCheck::default();
        let beta = F::random(OsRng);
        let gamma = F::random(OsRng);
        let lookup_slices = lookups.iter().map(|lookup| lookup.as_slice()).collect_vec();
        let count_polys = lookup_check
            .commit_counts(&lookup_slices, &polys, &[], &beta)
            .unwrap();
        let h_polys = lookup_check.commit_hs(&count_polys, &gamma).unwrap();
        for h_poly in h_polys {
            assert_eq!(sum(h_poly.evals().to_vec().into_iter()), F::zero());
        }
        drop(lookup_slices);

        // Wrong path
        let mut lookup_check = LookupCheck::default();
        lookups[0][0] = (&poly_exps[0] + &poly_exps[1], &poly_exps[3] + &poly_exps[5]);
        let lookup_slices = lookups.iter().map(|lookup| lookup.as_slice()).collect_vec();
        assert!(lookup_check
            .commit_counts(&lookup_slices, &polys, &[], &beta)
            .is_err());
    }

    #[test]
    fn test_build_constraints() {
        let one = Expression::Constant(F::one());
        let two = Expression::Constant(F::one() + F::one());
        let lookups = vec![
            vec![(one.clone(), two.clone()), (two.clone(), one.clone())],
            vec![(one.clone(), one.clone()), (two.clone(), two.clone())],
        ];
        let circuit_info = {
            let mut circuit_info = PlonkishCircuitInfo::<F> {
                num_instances: vec![],
                num_witness_polys: vec![],
                num_challenges: vec![2],
                lookups,
                ..Default::default()
            };
            circuit_info.initialize_permutation_info(circuit_info.permutation_polys().as_slice());
            circuit_info
        };
        let count_offset = 0;
        let h_offset = 2;
        assert_eq!(circuit_info.lookup_count_offset(), count_offset);
        assert_eq!(circuit_info.lookup_h_offset(), h_offset);

        let beta = Expression::Constant(F::one() + F::one() + F::one());
        let gamma = Expression::Constant(F::one() + F::one() + F::one() + F::one());

        let counts: Vec<Expression<F>> = vec![
            Expression::Polynomial(Query::new(count_offset, Rotation::cur())),
            Expression::Polynomial(Query::new(count_offset + 1, Rotation::cur())),
        ];
        let hs: Vec<Expression<F>> = vec![
            Expression::Polynomial(Query::new(h_offset, Rotation::cur())),
            Expression::Polynomial(Query::new(h_offset + 1, Rotation::cur())),
        ];
        let lookups = &circuit_info.lookups;
        let fxs = vec![
            Expression::distribute_powers([&lookups[0][0].0, &lookups[0][1].0], &beta) + &gamma,
            Expression::distribute_powers([&lookups[1][0].0, &lookups[1][1].0], &beta) + &gamma,
        ];
        let gxs = vec![
            Expression::distribute_powers([&lookups[0][0].1, &lookups[0][1].1], &beta) + &gamma,
            Expression::distribute_powers([&lookups[1][0].1, &lookups[1][1].1], &beta) + &gamma,
        ];

        let (constraint_expressions, sumcheck_expressions) =
            LookupCheck::<F>::build_constraints(&circuit_info, &beta, &gamma);
        let want_constraint_expressions = (0..2)
            .map(|i| (&hs[i] * &gxs[i] + &counts[i]) * &fxs[i] - &gxs[i])
            .collect_vec();
        assert_eq!(constraint_expressions, want_constraint_expressions);
        assert_eq!(sumcheck_expressions, hs);
    }
}
