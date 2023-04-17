use halo2_curves::group::ff::PrimeField;

use crate::backend::{
    piop::sum_check::{
        classic::{ClassicSumCheck, EvaluationsProver},
        evaluate, lagrange_eval, SumCheck,
    },
    poly::multilinear::rotation_eval_points,
    util::{
        arithmetic::{inner_product, BooleanHypercube},
        expression::{Expression, Query, Rotation},
        transcript::FieldTranscriptWrite,
        Itertools,
    },
    Error,
};
use std::collections::{BTreeSet, HashMap};

#[allow(clippy::type_complexity)]
pub(super) fn verify_sum_check<F: PrimeField>(
    num_vars: usize,
    expression: &Expression<F>,
    instances: &[&[F]],
    sum_check_msgs: &[&[F]],
    challenges: &[F],
    y: &[F],
    sum: F,
    evals: &[F],
    transcript: &mut impl FieldTranscriptWrite<F>,
) -> Result<Vec<Vec<F>>, Error> {
    let (x_eval, x) = ClassicSumCheck::<EvaluationsProver<_, true>>::verify(
        &(),
        sum_check_msgs,
        num_vars,
        expression.degree(),
        sum,
        transcript,
    )?;

    let evals: HashMap<Query, F> = evals.iter().enumerate().map(|(i, e)| (Query::new(i, Rotation::cur()), *e)).collect();

    let pcs_query = pcs_query(expression, instances.len());
    let evals = instance_evals(num_vars, expression, instances, &x)
        .into_iter()
        .chain(evals)
        .collect();
    // println!("evals: {:?}", evals);
    if evaluate(expression, num_vars, &evals, challenges, &[y], &x) != x_eval {
        return Err(Error::InvalidSnark(
            "Unmatched between sum_check output and query evaluation".to_string(),
        ));
    }
    Ok(points(&pcs_query, &x))
}

fn instance_evals<F: PrimeField>(
    num_vars: usize,
    expression: &Expression<F>,
    instances: &[&[F]],
    x: &[F],
) -> Vec<(Query, F)> {
    let mut instance_query = expression.used_query();
    instance_query.retain(|query| query.poly() < instances.len());

    let lagranges = {
        let mut lagranges = instance_query.iter().fold(0..0, |range, query| {
            let i = -query.rotation().0;
            range.start.min(i)..range.end.max(i + instances[query.poly()].len() as i32)
        });
        if lagranges.start < 0 {
            lagranges.start -= 1;
        }
        if lagranges.end > 0 {
            lagranges.end += 1;
        }
        lagranges
    };

    let bh = BooleanHypercube::new(num_vars).iter().collect_vec();
    let lagrange_evals = lagranges
        .filter_map(|i| {
            (i != 0).then(|| {
                let b = bh[i.rem_euclid(1 << num_vars as i32) as usize];
                (i, lagrange_eval(x, b))
            })
        })
        .collect::<HashMap<_, _>>();

    instance_query
        .into_iter()
        .map(|query| {
            let is = if query.rotation() > Rotation::cur() {
                (-query.rotation().0..0)
                    .chain(1..)
                    .take(instances[query.poly()].len())
                    .collect_vec()
            } else {
                (1 - query.rotation().0..)
                    .take(instances[query.poly()].len())
                    .collect_vec()
            };
            let eval = inner_product(
                instances[query.poly()],
                is.iter().map(|i| lagrange_evals.get(i).unwrap()),
            );
            (query, eval)
        })
        .collect()
}

pub(super) fn pcs_query<F: PrimeField>(
    expression: &Expression<F>,
    num_instance_poly: usize,
) -> BTreeSet<Query> {
    let mut used_query = expression.used_query();
    used_query.retain(|query| query.poly() >= num_instance_poly);
    used_query
}

pub(super) fn points<F: PrimeField>(pcs_query: &BTreeSet<Query>, x: &[F]) -> Vec<Vec<F>> {
    pcs_query
        .iter()
        .map(Query::rotation)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .flat_map(|rotation| rotation_eval_points(x, rotation))
        .collect_vec()
}

pub(super) fn point_offset(pcs_query: &BTreeSet<Query>) -> HashMap<Rotation, usize> {
    let rotations = pcs_query
        .iter()
        .map(Query::rotation)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect_vec();
    rotations.windows(2).fold(
        HashMap::from_iter([(rotations[0], 0)]),
        |mut point_offset, rotations| {
            let last_rotation = rotations[0];
            let offset = point_offset[&last_rotation] + (1 << last_rotation.distance());
            point_offset.insert(rotations[1], offset);
            point_offset
        },
    )
}
