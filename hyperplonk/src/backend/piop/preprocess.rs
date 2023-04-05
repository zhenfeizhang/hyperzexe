use halo2_curves::group::ff::PrimeField;

use crate::backend::{
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::steps,
        expression::{Expression, Query, Rotation},
        Itertools,
    },
    PlonkishCircuitInfo,
};
use core::hash::Hash;
use std::{array, iter, mem, sync::Arc};

use super::{lookup_check::LookupCheck, perm_check::PermutationCheck};

pub(super) fn generate_preprocessed_poly_group<F: PrimeField>(
    circuit_info: &PlonkishCircuitInfo<F>,
) -> (
    Vec<Arc<MultilinearPolynomial<F>>>,
    Vec<Arc<MultilinearPolynomial<F>>>,
) {
    assert!(circuit_info.permute_info_initialized);
    let num_vars = circuit_info.k;
    // Compute preprocesses
    let preprocess_polys = circuit_info
        .preprocess_polys
        .iter()
        .cloned()
        .map(|p| Arc::new(MultilinearPolynomial::new(p)))
        .collect_vec();

    // Compute permutation polys
    let permutation_polys = permutation_polys(
        num_vars,
        &circuit_info.permutation_polys(),
        &circuit_info.permutations,
    );
    (preprocess_polys, permutation_polys)
}

// The order of polynomials are
// instance polynomials
// preprocessed polynomials
// witness polynomials
// permutation polynomials
// lookup count polynomials
// lookup h polynomials
// permutation frac polynomials
// permutation product polynomials
// permutation p1 and p2 polynomials
pub(super) fn compose<F>(
    circuit_info: &PlonkishCircuitInfo<F>,
) -> (usize, Expression<F>, Expression<F>)
where
    F: PrimeField + Hash,
{
    debug_assert!(circuit_info.permute_info_initialized);
    let challenge_offset = circuit_info.num_challenges.iter().sum::<usize>();
    let [beta, gamma, alpha, eta] =
        &array::from_fn(|idx| Expression::<F>::Challenge(challenge_offset + idx));
    let (lookup_constraints, lookup_sumchecks) =
        LookupCheck::build_constraints(circuit_info, beta, gamma);

    let max_degree = iter::empty()
        .chain(circuit_info.constraints.iter())
        .chain(lookup_constraints.iter())
        .map(Expression::degree)
        .chain(circuit_info.max_degree)
        .chain(Some(2))
        .max()
        .unwrap();
    let permutation_constraints = PermutationCheck::build_constraints(circuit_info, beta, gamma);

    let constraint_expression = {
        let constraints = circuit_info
            .constraints
            .iter()
            .chain(lookup_constraints.iter())
            .chain(permutation_constraints.iter())
            .collect_vec();
        let eq = Expression::eq_xy(0);
        let c1 = Expression::distribute_powers(constraints, alpha) * eq;
        Expression::distribute_powers(iter::once(&c1).chain(lookup_sumchecks.iter()), alpha)
    };
    let frac_left_eq = &Expression::<F>::eq_xy_shifted_00(0);
    let frac_right_eq = &Expression::<F>::eq_xy_shifted_01(0);
    let prod_left_eq = &Expression::<F>::eq_xy_shifted_10(0);
    let prod_right_eq = &Expression::<F>::eq_xy_shifted_11(0);

    let num_chunks = circuit_info.num_permutation_chunks();
    let opening_expression: Expression<F> = {
        let other_polys = (circuit_info.preprocess_offset()
            ..circuit_info.permutation_p1_p2_offset())
            .map(|poly| {
                let poly = Query::new(poly, Rotation::cur());
                Expression::<F>::Polynomial(poly)
            })
            .collect_vec();
        let permute_frac_polys = (circuit_info.permutation_frac_offset()
            ..circuit_info.permutation_frac_offset() + num_chunks)
            .map(|poly| {
                let poly = Query::new(poly, Rotation::cur());
                Expression::Polynomial(poly)
            });
        let permute_prod_polys = (circuit_info.permutation_prod_offset()
            ..circuit_info.permutation_prod_offset() + num_chunks)
            .map(|poly| {
                let poly = Query::new(poly, Rotation::cur());
                Expression::Polynomial(poly)
            });
        let p1_p2_polys = permute_frac_polys
            .zip(permute_prod_polys)
            .flat_map(|(frac, prod)| {
                let p1 = frac_left_eq * &frac + prod_left_eq * &prod;
                let p2 = frac_right_eq * &frac + prod_right_eq * &prod;
                [p1, p2]
            });
        let single_point_eq: Expression<F> = Expression::eq_xy(0);
        let opening_expression_normal_point =
            Expression::distribute_powers(other_polys.as_slice(), eta);
        let qs = iter::once(opening_expression_normal_point * &single_point_eq)
            .chain(p1_p2_polys)
            .collect_vec();
        Expression::distribute_powers(qs.iter(), eta)
    };

    (max_degree, constraint_expression, opening_expression)
}

pub(super) fn permutation_polys<F: PrimeField>(
    num_vars: usize,
    permutation_polys: &[usize],
    cycles: &[Vec<(usize, usize)>],
) -> Vec<Arc<MultilinearPolynomial<F>>> {
    let poly_index = {
        let mut poly_index = vec![0; permutation_polys.last().map(|poly| 1 + poly).unwrap_or(0)];
        for (idx, poly) in permutation_polys.iter().enumerate() {
            poly_index[*poly] = idx;
        }
        poly_index
    };
    let mut permutations = (0..permutation_polys.len() as u64)
        .map(|idx| {
            steps(F::from(idx << num_vars))
                .take(1 << num_vars)
                .collect_vec()
        })
        .collect_vec();
    for cycle in cycles.iter() {
        let (i0, j0) = cycle[0];
        let mut last = permutations[poly_index[i0]][j0];
        for &(i, j) in cycle.iter().cycle().skip(1).take(cycle.len()) {
            assert_ne!(j, 0);
            mem::swap(&mut permutations[poly_index[i]][j], &mut last);
        }
    }
    permutations
        .into_iter()
        .map(|permutation| Arc::new(MultilinearPolynomial::new(permutation)))
        .collect()
}

#[cfg(test)]
pub(crate) mod test {
    use crate::backend::{
        piop::util::{plonk_expression, plonk_with_lookup_expression},
        util::expression::{Expression, Query, Rotation},
    };
    use halo2_curves::bn256::Fr;
    use std::array;

    #[test]
    fn compose_plonk() {
        let (constraint_expression, opening_expression) = plonk_expression();
        assert_eq!((constraint_expression, opening_expression), {
            let [pi, q_l, q_r, q_m, q_o, q_c, w_l, w_r, w_o, s_1, s_2, s_3] =
                &array::from_fn(|poly| Query::new(poly, Rotation::cur()))
                    .map(Expression::Polynomial);
            let [frac, prod] = &[
                Query::new(12, Rotation::cur()),
                Query::new(13, Rotation::cur()),
            ]
            .map(Expression::Polynomial);
            let [p1, p2] = &[
                Query::new(14, Rotation::cur()),
                Query::new(15, Rotation::cur()),
            ]
            .map(Expression::Polynomial);
            let [beta, gamma, alpha, eta] = &[0, 1, 2, 3].map(Expression::<Fr>::Challenge);
            let [id_1, id_2, id_3] = array::from_fn(Expression::identity);
            let l_minus_2 = Expression::<Fr>::lagrange(-2);
            let one = Expression::Constant(Fr::one());
            let constraints = {
                vec![
                    q_l * w_l + q_r * w_r + q_m * w_l * w_r + q_o * w_o + q_c + pi,
                    prod - p1 * p2,
                    frac * ((w_l + beta * s_1 + gamma)
                        * (w_r + beta * s_2 + gamma)
                        * (w_o + beta * s_3 + gamma))
                        - ((w_l + beta * id_1 + gamma)
                            * (w_r + beta * id_2 + gamma)
                            * (w_o + beta * id_3 + gamma)),
                    l_minus_2 * (prod - one),
                ]
            };
            let eq = &Expression::eq_xy(0);
            let frac_left_eq: Expression<Fr> = Expression::eq_xy_shifted_00(0);
            let frac_right_eq: Expression<Fr> = Expression::eq_xy_shifted_01(0);
            let prod_left_eq: Expression<Fr> = Expression::eq_xy_shifted_10(0);
            let prod_right_eq: Expression<Fr> = Expression::eq_xy_shifted_11(0);
            let constraint_expression = Expression::distribute_powers(
                &[Expression::distribute_powers(&constraints, &alpha) * eq],
                &alpha,
            );
            let open_p1 = frac_left_eq * frac + prod_left_eq * prod;
            let open_p2 = frac_right_eq * frac + prod_right_eq * prod;
            let openings = vec![
                q_l, q_r, q_m, q_o, q_c, w_l, w_r, w_o, s_1, s_2, s_3, frac, prod,
            ];
            let opening_expression_normal_point = Expression::distribute_powers(openings, eta) * eq;
            let opening_expression = Expression::distribute_powers(
                &[opening_expression_normal_point, open_p1, open_p2],
                eta,
            );
            (constraint_expression, opening_expression)
        });
    }

    #[test]
    fn compose_plonk_with_lookup() {
        let (constraint_expression, opening_expression) = plonk_with_lookup_expression();
        assert_eq!((constraint_expression, opening_expression), {
            let [pi, q_l, q_r, q_m, q_o, q_c, q_lookup, t_l, t_r, t_o, w_l, w_r, w_o, s_1, s_2, s_3] =
                &array::from_fn(|poly| Query::new(poly, Rotation::cur()))
                    .map(Expression::Polynomial);
            let [lookup_count, lookup_h] = &[
                Query::new(16, Rotation::cur()),
                Query::new(17, Rotation::cur()),
            ]
            .map(Expression::Polynomial);
            let [perm_frac, perm_prod] = &[
                Query::new(18, Rotation::cur()),
                Query::new(19, Rotation::cur()),
            ]
            .map(Expression::Polynomial);
            let [p1, p2] = &[
                Query::new(20, Rotation::cur()),
                Query::new(21, Rotation::cur()),
            ]
            .map(Expression::Polynomial);
            let [beta, gamma, alpha, eta] = &array::from_fn(Expression::<Fr>::Challenge);
            let [id_1, id_2, id_3] = array::from_fn(Expression::identity);
            let l_minus_2 = Expression::<Fr>::lagrange(-2);
            let one = &Expression::Constant(Fr::one());
            let lookup_compressed_input =
                Expression::distribute_powers(&[w_l, w_r, w_o].map(|w| q_lookup * w), beta) + gamma;
            let lookup_compressed_table =
                Expression::distribute_powers([t_l, t_r, t_o], beta) + gamma;
            let constraints = {
                vec![
                    q_l * w_l + q_r * w_r + q_m * w_l * w_r + q_o * w_o + q_c + pi,
                    (lookup_h * &lookup_compressed_table + lookup_count) * lookup_compressed_input
                        - &lookup_compressed_table,
                    perm_prod - p1 * p2,
                    perm_frac
                        * ((w_l + beta * s_1 + gamma)
                            * (w_r + beta * s_2 + gamma)
                            * (w_o + beta * s_3 + gamma))
                        - ((w_l + beta * id_1 + gamma)
                            * (w_r + beta * id_2 + gamma)
                            * (w_o + beta * id_3 + gamma)),
                    l_minus_2 * (perm_prod - one),
                ]
            };
            let eq = &Expression::eq_xy(0);
            let frac_left_eq: Expression<Fr> = Expression::eq_xy_shifted_00(0);
            let frac_right_eq: Expression<Fr> = Expression::eq_xy_shifted_01(0);
            let prod_left_eq: Expression<Fr> = Expression::eq_xy_shifted_10(0);
            let prod_right_eq: Expression<Fr> = Expression::eq_xy_shifted_11(0);
            let constraint_expression = Expression::distribute_powers(
                &[
                    Expression::distribute_powers(&constraints, &alpha) * eq,
                    lookup_h.clone(),
                ],
                &alpha,
            );
            let open_p1 = frac_left_eq * perm_frac + prod_left_eq * perm_prod;
            let open_p2 = frac_right_eq * perm_frac + prod_right_eq * perm_prod;
            let openings = vec![
                q_l,
                q_r,
                q_m,
                q_o,
                q_c,
                q_lookup,
                t_l,
                t_r,
                t_o,
                w_l,
                w_r,
                w_o,
                s_1,
                s_2,
                s_3,
                lookup_count,
                lookup_h,
                perm_frac,
                perm_prod,
            ];
            let opening_expression_normal_point = Expression::distribute_powers(openings, eta) * eq;
            let opening_expression = Expression::distribute_powers(
                &[opening_expression_normal_point, open_p1, open_p2],
                eta,
            );
            (constraint_expression, opening_expression)
        });
    }
}
