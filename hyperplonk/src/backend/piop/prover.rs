use halo2_curves::group::ff::Field;

use crate::backend::{
    pcs::Evaluation,
    piop::{
        sum_check::{
            classic::{ClassicSumCheck, EvaluationsProver},
            SumCheck, VirtualPolynomial,
        },
        verifier::{pcs_query, point_offset, points},
    },
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::{BooleanHypercube, PrimeField},
        end_timer,
        expression::{Expression, Rotation},
        start_timer,
        transcript::FieldTranscriptWrite,
        Itertools,
    },
    Error,
};
use std::{iter::Chain, slice::Iter, sync::Arc};

pub(super) fn instances_polys<'a, F: PrimeField>(
    num_vars: usize,
    instances: impl IntoIterator<Item = impl IntoIterator<Item = &'a F>>,
) -> Vec<Arc<MultilinearPolynomial<F>>> {
    let bh = BooleanHypercube::new(num_vars);
    instances
        .into_iter()
        .map(|instances| {
            let mut poly = vec![F::zero(); 1 << num_vars];
            for (b, instance) in bh.iter().skip(1).zip(instances.into_iter()) {
                poly[b] = *instance;
            }
            poly
        })
        .map(|p| Arc::new(MultilinearPolynomial::new(p)))
        .collect()
}

#[allow(clippy::type_complexity)]
pub(super) fn prove_sum_check<F: PrimeField>(
    num_instance_poly: usize,
    expression: &Expression<F>,
    polys: &[Arc<MultilinearPolynomial<F>>],
    challenges: &[F],
    y: Vec<F>,
    sum: F,
    transcript: &mut impl FieldTranscriptWrite<F>,
) -> Result<(Vec<Vec<F>>, Vec<Vec<F>>, Vec<Evaluation<F>>), Error> {
    if polys.is_empty() {
        return Ok((vec![], vec![], vec![]));
    }
    let num_vars = polys[0].num_vars();
    let ys = [y];

    let virtual_poly = VirtualPolynomial::new(expression, polys.to_vec(), &challenges, &ys);
    let (sum_check_proof, x, evals) = ClassicSumCheck::<EvaluationsProver<_, true>>::prove(
        &(),
        num_vars,
        virtual_poly,
        sum,
        transcript,
    )?;

    let pcs_query = pcs_query(expression, num_instance_poly);
    let point_offset = point_offset(&pcs_query);

    let timer = start_timer(|| format!("evals-{}", pcs_query.len()));
    let evals = pcs_query
        .iter()
        .flat_map(|query| {
            (point_offset[&query.rotation()]..)
                .zip(if query.rotation() == Rotation::cur() {
                    vec![evals[query.poly()]]
                } else {
                    polys[query.poly()].evaluate_for_rotation(&x, query.rotation())
                })
                .map(|(point, eval)| Evaluation::new(query.poly(), point, eval))
        })
        .collect_vec();
    end_timer(timer);

    transcript.write_field_elements(evals.iter().map(Evaluation::value))?;
    Ok((sum_check_proof, points(&pcs_query, &x), evals))
}

pub(super) fn extend_polys_and_compute_offsets<F>(
    polys: &mut Vec<Arc<MultilinearPolynomial<F>>>,
    permutation: &[Arc<MultilinearPolynomial<F>>],
    lookup_count: &[Arc<MultilinearPolynomial<F>>],
    lookup_h: &[Arc<MultilinearPolynomial<F>>],
    permutation_frac: &[Arc<MultilinearPolynomial<F>>],
    permutation_prod: &[Arc<MultilinearPolynomial<F>>],
    p1_p2: &[Arc<MultilinearPolynomial<F>>],
) -> (usize, usize, usize, usize, usize, usize) {
    let mut extend_polys = |others: Vec<Arc<MultilinearPolynomial<F>>>| {
        let offset = polys.len();
        polys.extend(others);
        offset
    };
    let permutation_offset = extend_polys(permutation.to_vec());
    let lookup_count_offset = extend_polys(lookup_count.to_vec());
    let lookup_h_offset = extend_polys(lookup_h.to_vec());
    let permute_frac_offset = extend_polys(permutation_frac.to_vec());
    let permute_prod_offset = extend_polys(permutation_prod.to_vec());
    let p1_p2_offset = extend_polys(p1_p2.to_vec());
    (
        permutation_offset,
        lookup_count_offset,
        lookup_h_offset,
        permute_frac_offset,
        permute_prod_offset,
        p1_p2_offset,
    )
}

pub(super) fn reorder_into_groups<T: Clone>(
    arr: Vec<T>,
    segment_groups: &[Vec<(usize, usize)>], // offset and length
) -> Vec<Vec<T>> {
    let mut groups = vec![];
    for segments in segment_groups {
        let mut group = vec![];
        for (offset, length) in segments {
            group.extend(arr[*offset..*offset + *length].to_vec());
        }
        groups.push(group);
    }
    groups
}

pub trait ArcItertools<T, U> {
    fn collect_vec_arc(self) -> Vec<Arc<T>>;
}

impl<I, F> ArcItertools<MultilinearPolynomial<F>, MultilinearPolynomial<F>> for I
where
    F: Field,
    I: Iterator<Item = MultilinearPolynomial<F>>,
{
    fn collect_vec_arc(self) -> Vec<Arc<MultilinearPolynomial<F>>>
    where
        I: Iterator<Item = MultilinearPolynomial<F>>,
    {
        self.map(Arc::new).collect_vec()
    }
}

impl<I, U, F> ArcItertools<MultilinearPolynomial<F>, U> for I
where
    F: Field,
    I: Iterator<Item = U>,
    U: Iterator<Item = MultilinearPolynomial<F>>,
{
    fn collect_vec_arc(self) -> Vec<Arc<MultilinearPolynomial<F>>> {
        self.flat_map(|p| p.into_iter().map(Arc::new)).collect_vec()
    }
}

pub trait ArcItertoolsDeref<T> {
    fn collect_vec_arc(self) -> Vec<Arc<T>>;
}

impl<F> ArcItertoolsDeref<MultilinearPolynomial<F>> for Iter<'_, MultilinearPolynomial<F>>
where
    F: Field,
{
    fn collect_vec_arc(self) -> Vec<Arc<MultilinearPolynomial<F>>> {
        self.map(|poly| Arc::new(poly.clone())).collect_vec()
    }
}

impl<F> ArcItertoolsDeref<MultilinearPolynomial<F>>
    for Chain<Iter<'_, MultilinearPolynomial<F>>, Iter<'_, MultilinearPolynomial<F>>>
where
    F: Field,
{
    fn collect_vec_arc(self) -> Vec<Arc<MultilinearPolynomial<F>>> {
        self.map(|poly| Arc::new(poly.clone())).collect_vec()
    }
}
