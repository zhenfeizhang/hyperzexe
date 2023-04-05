use crate::{
    backend::util::{
        arithmetic::{div_ceil, usize_from_bits_be, BooleanHypercube, Field},
        expression::Rotation,
        impl_index,
        parallel::{num_threads, parallelize, parallelize_iter},
        BitIndex, Itertools,
    },
    Error,
};
use num_integer::Integer;
use rand::RngCore;
use rayon::prelude::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use std::{
    borrow::Cow,
    iter::{self, Sum},
    mem,
    ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign},
    sync::Arc,
};

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct MultilinearPolynomial<F> {
    evals: Vec<F>,
    num_vars: usize,
}

impl<F: Field> Default for MultilinearPolynomial<F> {
    fn default() -> Self {
        MultilinearPolynomial::zero()
    }
}

impl<F: Field> MultilinearPolynomial<F> {
    pub fn new(evals: Vec<F>) -> Self {
        let num_vars = if evals.is_empty() {
            0
        } else {
            let num_vars = evals.len().ilog2() as usize;
            debug_assert_eq!(evals.len(), 1 << num_vars);
            num_vars
        };

        Self { evals, num_vars }
    }

    pub fn zero() -> Self {
        Self {
            evals: vec![F::one()],
            num_vars: 0,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.num_vars == 0
    }

    pub fn evals(&self) -> &[F] {
        self.evals.as_slice()
    }

    pub fn into_evals(self) -> Vec<F> {
        self.evals
    }

    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn iter(&self) -> impl Iterator<Item = &F> {
        self.evals.iter()
    }
}

impl<F: Field> MultilinearPolynomial<F> {
    pub fn eq_xy(y: &[F]) -> Self {
        if y.is_empty() {
            return Self::zero();
        }

        let expand_serial = |next_evals: &mut [F], evals: &[F], y_i: &F| {
            for (next_evals, eval) in next_evals.chunks_mut(2).zip(evals.iter()) {
                next_evals[1] = *eval * y_i;
                next_evals[0] = *eval - &next_evals[1];
            }
        };

        let mut evals = vec![F::one()];
        for y_i in y.iter().rev() {
            let mut next_evals = vec![F::zero(); 2 * evals.len()];
            if evals.len() < 32 {
                expand_serial(&mut next_evals, &evals, y_i);
            } else {
                let mut chunk_size = div_ceil(evals.len(), num_threads());
                if chunk_size.is_odd() {
                    chunk_size += 1;
                }
                parallelize_iter(
                    next_evals
                        .chunks_mut(chunk_size)
                        .zip(evals.chunks(chunk_size >> 1)),
                    |(next_evals, evals)| expand_serial(next_evals, evals, y_i),
                );
            }
            evals = next_evals;
        }

        Self {
            evals,
            num_vars: y.len(),
        }
    }

    /// Compute ((1 - yn) * (1 - a) + yn * a ) eq ((b, y1, ..., yn), x)
    pub fn eq_xy_shifted_ab(y: &[F], a: bool, b: bool) -> Self {
        let y_last = if a {
            y[y.len() - 1]
        } else {
            F::one() - y[y.len() - 1]
        };
        let y_0 = if b { F::one() } else { F::zero() };
        let y = iter::once(y_0)
            .chain(y[..y.len() - 1].to_vec())
            .collect_vec();
        let mut res = Self::eq_xy(&y).into_evals();
        res.iter_mut().for_each(|eval| *eval *= &y_last);
        Self::new(res)
    }

    pub fn rand(num_vars: usize, mut rng: impl RngCore) -> Self {
        Self::new(
            iter::repeat_with(|| F::random(&mut rng))
                .take(1 << num_vars)
                .collect(),
        )
    }

    pub fn evaluate(&self, x: &[F]) -> F {
        debug_assert_eq!(x.len(), self.num_vars);

        let mut evals = Cow::Borrowed(self.evals());
        let mut bits = Vec::new();
        let mut buf = Vec::with_capacity(self.evals.len() >> 1);
        for x_i in x.iter() {
            if x_i == &F::zero() || x_i == &F::one() {
                bits.push(x_i == &F::one());
                continue;
            }

            let distance = bits.len() + 1;
            let skip = usize_from_bits_be(bits.drain(..).rev());
            merge(&mut evals, x_i, distance, skip, &mut buf);
        }

        evals[usize_from_bits_be(bits.drain(..).rev())]
    }

    pub fn fix_variable_into(&self, target: &mut Self, x_i: &F) {
        target.num_vars = self.num_vars - 1;
        merge_into(self.evals(), &mut target.evals, x_i, 1, 0);
    }

    pub fn fix_variable(&mut self, x_i: &F, buf: &mut Self) {
        self.fix_variable_into(buf, x_i);
        mem::swap(self, buf);
    }

    pub fn fix_variables(&self, partial_point: &[F]) -> MultilinearPolynomial<F> {
        debug_assert!(
            partial_point.len() <= self.num_vars,
            "invalid size of partial point"
        );
        let nv = self.num_vars;
        let mut poly = self.clone();
        let mut buf = MultilinearPolynomial::new(vec![F::zero(); 1 << (nv - 1)]);
        for (i, x_i) in partial_point.iter().enumerate() {
            buf.num_vars = nv - i - 1;
            merge_into(poly.evals(), &mut buf.evals, x_i, 1, 0);
            mem::swap(&mut poly, &mut buf);
        }
        return poly;
    }

    /// merge a set of polynomials. Returns an error if the
    /// polynomials do not share a same number of nvs.
    pub fn merge_polynomials(polynomials: &[Arc<Self>]) -> Result<Arc<Self>, Error> {
        let mut scalars = vec![];
        for poly in polynomials.iter() {
            scalars.extend_from_slice(poly.evals());
        }
        scalars.extend_from_slice(
            vec![F::zero(); (scalars.len().next_power_of_two()) - scalars.len()].as_ref(),
        );
        Ok(Arc::new(Self::new(scalars)))
    }

    pub fn fix_last_variables(&self, partial_point: &[F]) -> MultilinearPolynomial<F> {
        debug_assert!(
            partial_point.len() <= self.num_vars,
            "invalid size of partial point"
        );
        let nv = self.num_vars;
        let mut poly = self.evals().to_vec();
        let dim = partial_point.len();
        // evaluate single variable of partial point from left to right
        for (i, point) in partial_point.iter().rev().enumerate().take(dim) {
            poly = fix_last_variable_helper(&poly, nv - i, point);
        }

        MultilinearPolynomial::<F>::new(poly)
    }

    pub fn evaluate_for_rotation(&self, x: &[F], rotation: Rotation) -> Vec<F> {
        debug_assert_eq!(x.len(), self.num_vars);
        if rotation == Rotation::cur() {
            return vec![self.evaluate(x)];
        }

        let distance = rotation.distance();
        let num_x = self.num_vars - distance;
        let mut evals = vec![F::zero(); 1 << distance];
        let chunk_size = div_ceil(evals.len(), num_threads());
        if rotation < Rotation::cur() {
            let x = &x[distance..];
            let flipped_x = x.iter().map(flip).collect_vec();
            let pattern = rotation_eval_point_pattern::<false>(self.num_vars, distance);
            let offset_mask = (1 << self.num_vars) - (1 << num_x);
            parallelize_iter(
                evals.chunks_mut(chunk_size).zip(pattern.chunks(chunk_size)),
                |(evals, pattern)| {
                    let mut buf = vec![F::zero(); 1 << (num_x - 1)];
                    let mut last_buf = Some(vec![F::zero(); 1 << (num_x - 1)]);
                    for (eval, pat) in evals.iter_mut().zip(pattern.iter()) {
                        let offset = pat & offset_mask;
                        let mut evals = Cow::Borrowed(&self[offset..offset + (1 << num_x)]);
                        for (idx, (x_i, flipped_x_i)) in x.iter().zip(flipped_x.iter()).enumerate()
                        {
                            let x_i = if pat.nth_bit(idx) { flipped_x_i } else { x_i };
                            merge_into(&evals, &mut buf, x_i, 1, 0);
                            if let Cow::Owned(_) = evals {
                                mem::swap(evals.to_mut(), &mut buf);
                            } else {
                                evals = mem::replace(&mut buf, last_buf.take().unwrap()).into();
                            }
                        }
                        *eval = evals[0];
                        last_buf = Some(evals.into_owned());
                    }
                },
            );
        } else {
            let x = &x[..num_x];
            let flipped_x = x.iter().map(flip).collect_vec();
            let pattern = rotation_eval_point_pattern::<true>(self.num_vars, distance);
            let skip_mask = (1 << distance) - 1;
            parallelize_iter(
                evals.chunks_mut(chunk_size).zip(pattern.chunks(chunk_size)),
                |(evals, pattern)| {
                    let mut buf = Vec::with_capacity(self.evals.len() >> 1);
                    let mut last_buf = Some(Vec::with_capacity(self.evals.len() >> 1));
                    for (eval, pat) in evals.iter_mut().zip(pattern.iter()) {
                        let mut evals = Cow::Borrowed(self.evals());
                        let skip = pat & skip_mask;
                        let x_0 = if pat.nth_bit(distance) {
                            &flipped_x[0]
                        } else {
                            &x[0]
                        };
                        merge_into(&evals, &mut buf, x_0, distance + 1, skip);
                        evals = mem::replace(&mut buf, last_buf.take().unwrap()).into();

                        for ((x_i, flipped_x_i), idx) in
                            x.iter().zip(flipped_x.iter()).zip(distance..).skip(1)
                        {
                            let x_i = if pat.nth_bit(idx) { flipped_x_i } else { x_i };
                            merge(&mut evals, x_i, 1, 0, &mut buf);
                        }
                        *eval = evals[0];
                        last_buf = Some(evals.into_owned());
                    }
                },
            );
        }
        evals
    }
}

impl<'lhs, 'rhs, F: Field> Add<&'rhs MultilinearPolynomial<F>> for &'lhs MultilinearPolynomial<F> {
    type Output = MultilinearPolynomial<F>;

    fn add(self, rhs: &'rhs MultilinearPolynomial<F>) -> MultilinearPolynomial<F> {
        let mut output = self.clone();
        output += rhs;
        output
    }
}

impl<'rhs, F: Field> AddAssign<&'rhs MultilinearPolynomial<F>> for MultilinearPolynomial<F> {
    fn add_assign(&mut self, rhs: &'rhs MultilinearPolynomial<F>) {
        debug_assert_eq!(self.num_vars, rhs.num_vars);

        parallelize(&mut self.evals, |(lhs, start)| {
            for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                *lhs += rhs;
            }
        });
    }
}

impl<'rhs, F: Field> AddAssign<(&'rhs F, &'rhs MultilinearPolynomial<F>)>
    for MultilinearPolynomial<F>
{
    fn add_assign(&mut self, (scalar, rhs): (&'rhs F, &'rhs MultilinearPolynomial<F>)) {
        debug_assert_eq!(self.num_vars, rhs.num_vars);

        if scalar == &F::one() {
            *self += rhs;
        } else if scalar != &F::zero() {
            parallelize(&mut self.evals, |(lhs, start)| {
                for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                    *lhs += &(*scalar * rhs);
                }
            });
        }
    }
}

impl<'lhs, 'rhs, F: Field> Sub<&'rhs MultilinearPolynomial<F>> for &'lhs MultilinearPolynomial<F> {
    type Output = MultilinearPolynomial<F>;

    fn sub(self, rhs: &'rhs MultilinearPolynomial<F>) -> MultilinearPolynomial<F> {
        let mut output = self.clone();
        output -= rhs;
        output
    }
}

impl<'rhs, F: Field> SubAssign<&'rhs MultilinearPolynomial<F>> for MultilinearPolynomial<F> {
    fn sub_assign(&mut self, rhs: &'rhs MultilinearPolynomial<F>) {
        debug_assert_eq!(self.num_vars, rhs.num_vars);

        parallelize(&mut self.evals, |(lhs, start)| {
            for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                *lhs -= rhs;
            }
        });
    }
}

impl<'lhs, 'rhs, F: Field> Mul<&'rhs F> for &'lhs MultilinearPolynomial<F> {
    type Output = MultilinearPolynomial<F>;

    fn mul(self, rhs: &'rhs F) -> MultilinearPolynomial<F> {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}

impl<'rhs, F: Field> MulAssign<&'rhs F> for MultilinearPolynomial<F> {
    fn mul_assign(&mut self, rhs: &'rhs F) {
        if rhs == &F::zero() {
            self.evals = vec![F::zero(); self.evals.len()]
        } else if rhs != &F::one() {
            parallelize(&mut self.evals, |(lhs, _)| {
                for lhs in lhs.iter_mut() {
                    *lhs *= rhs;
                }
            });
        }
    }
}

impl<'a, F: Field> Sum<&'a MultilinearPolynomial<F>> for MultilinearPolynomial<F> {
    fn sum<I: Iterator<Item = &'a MultilinearPolynomial<F>>>(
        mut iter: I,
    ) -> MultilinearPolynomial<F> {
        let init = match (iter.next(), iter.next()) {
            (Some(lhs), Some(rhs)) => lhs + rhs,
            (Some(lhs), None) => return lhs.clone(),
            _ => unreachable!(),
        };
        iter.fold(init, |mut acc, poly| {
            acc += poly;
            acc
        })
    }
}

impl<F: Field> Sum<MultilinearPolynomial<F>> for MultilinearPolynomial<F> {
    fn sum<I: Iterator<Item = MultilinearPolynomial<F>>>(iter: I) -> MultilinearPolynomial<F> {
        iter.reduce(|mut acc, poly| {
            acc += &poly;
            acc
        })
        .unwrap()
    }
}

impl<F: Field> Sum<(F, MultilinearPolynomial<F>)> for MultilinearPolynomial<F> {
    fn sum<I: Iterator<Item = (F, MultilinearPolynomial<F>)>>(
        mut iter: I,
    ) -> MultilinearPolynomial<F> {
        let (scalar, mut poly) = iter.next().unwrap();
        poly *= &scalar;
        iter.fold(poly, |mut acc, (scalar, poly)| {
            acc += (&scalar, &poly);
            acc
        })
    }
}

impl_index!(MultilinearPolynomial, evals);

pub fn rotation_eval<F: Field>(x: &[F], rotation: Rotation, evals_for_rotation: &[F]) -> F {
    if rotation == Rotation::cur() {
        assert!(evals_for_rotation.len() == 1);
        return evals_for_rotation[0];
    }

    let num_vars = x.len();
    let distance = rotation.distance();
    assert!(evals_for_rotation.len() == 1 << distance);
    assert!(distance <= num_vars);

    let (pattern, nths, x) = if rotation < Rotation::cur() {
        (
            rotation_eval_coeff_pattern::<false>(num_vars, distance),
            (1..=distance).rev().collect_vec(),
            x[0..distance].iter().rev().collect_vec(),
        )
    } else {
        (
            rotation_eval_coeff_pattern::<true>(num_vars, distance),
            (num_vars - 1..).take(distance).collect(),
            x[num_vars - distance..].iter().collect(),
        )
    };
    x.into_iter().zip(nths).enumerate().fold(
        Cow::Borrowed(evals_for_rotation),
        |evals, (idx, (x_i, nth))| {
            pattern
                .iter()
                .step_by(1 << idx)
                .map(|pat| pat.nth_bit(nth))
                .zip(zip_self!(evals.iter()))
                .map(|(bit, (eval_0, eval_1))| {
                    if bit {
                        (*eval_0 - eval_1) * x_i + eval_1
                    } else {
                        (*eval_1 - eval_0) * x_i + eval_0
                    }
                })
                .collect_vec()
                .into()
        },
    )[0]
}

pub fn rotation_eval_points<F: Field>(x: &[F], rotation: Rotation) -> Vec<Vec<F>> {
    if rotation == Rotation::cur() {
        return vec![x.to_vec()];
    }

    let distance = rotation.distance();
    let num_x = x.len() - distance;
    if rotation < Rotation::cur() {
        let pattern = rotation_eval_point_pattern::<false>(x.len(), distance);
        let x = &x[distance..];
        let flipped_x = x.iter().map(flip).collect_vec();
        pattern
            .iter()
            .map(|pat| {
                iter::empty()
                    .chain((0..num_x).map(|idx| {
                        if pat.nth_bit(idx) {
                            flipped_x[idx]
                        } else {
                            x[idx]
                        }
                    }))
                    .chain((0..distance).map(|idx| bit_to_field(pat.nth_bit(idx + num_x))))
                    .collect_vec()
            })
            .collect()
    } else {
        let pattern = rotation_eval_point_pattern::<true>(x.len(), distance);
        let x = &x[..num_x];
        let flipped_x = x.iter().map(flip).collect_vec();
        pattern
            .iter()
            .map(|pat| {
                iter::empty()
                    .chain((0..distance).map(|idx| bit_to_field(pat.nth_bit(idx))))
                    .chain((0..num_x).map(|idx| {
                        if pat.nth_bit(idx + distance) {
                            flipped_x[idx]
                        } else {
                            x[idx]
                        }
                    }))
                    .collect_vec()
            })
            .collect()
    }
}

fn rotation_eval_point_pattern<const NEXT: bool>(num_vars: usize, distance: usize) -> Vec<usize> {
    let bh = BooleanHypercube::new(num_vars);
    let remainder = if NEXT { bh.primitive() } else { bh.x_inv() };
    let mut pattern = vec![0; 1 << distance];
    for depth in 0..distance {
        for (e, o) in zip_self!(0..pattern.len(), 1 << (distance - depth)) {
            let rotated = if NEXT {
                pattern[e] << 1
            } else {
                pattern[e] >> 1
            };
            pattern[o] = rotated ^ remainder;
            pattern[e] = rotated;
        }
    }
    pattern
}

fn rotation_eval_coeff_pattern<const NEXT: bool>(num_vars: usize, distance: usize) -> Vec<usize> {
    let bh = BooleanHypercube::new(num_vars);
    let remainder = if NEXT {
        bh.primitive() - (1 << num_vars)
    } else {
        bh.x_inv() << distance
    };
    let mut pattern = vec![0; 1 << (distance - 1)];
    for depth in 0..distance - 1 {
        for (e, o) in zip_self!(0..pattern.len(), 1 << (distance - depth - 1)) {
            let rotated = if NEXT {
                pattern[e] << 1
            } else {
                pattern[e] >> 1
            };
            pattern[o] = rotated ^ remainder;
            pattern[e] = rotated;
        }
    }
    pattern
}

fn flip<F: Field>(x: &F) -> F {
    F::one() - x
}

fn bit_to_field<F: Field>(bit: bool) -> F {
    if bit {
        F::one()
    } else {
        F::zero()
    }
}

fn merge<F: Field>(evals: &mut Cow<[F]>, x_i: &F, distance: usize, skip: usize, buf: &mut Vec<F>) {
    merge_into(evals, buf, x_i, distance, skip);
    if let Cow::Owned(_) = evals {
        mem::swap(evals.to_mut(), buf);
    } else {
        *evals = mem::replace(buf, Vec::with_capacity(buf.len() >> 1)).into();
    }
}

fn merge_into<F: Field>(evals: &[F], target: &mut Vec<F>, x_i: &F, distance: usize, skip: usize) {
    debug_assert!(target.capacity() >= evals.len() >> distance);
    target.resize_with(evals.len() >> distance, F::zero);

    let step = 1 << distance;
    parallelize(target, |(target, start)| {
        let start = (start << distance) + skip;
        for (target, (eval_0, eval_1)) in
            target.iter_mut().zip(zip_self!(evals.iter(), step, start))
        {
            *target = (*eval_1 - eval_0) * x_i + eval_0;
        }
    });
}

fn fix_last_variable_helper<F: Field>(data: &[F], nv: usize, point: &F) -> Vec<F> {
    let half_len = 1 << (nv - 1);
    let mut res = vec![F::zero(); half_len];

    // evaluate single variable of partial point from left to right
    #[cfg(not(feature = "parallel"))]
    for b in 0..half_len {
        res[b] = data[b] + (data[b + half_len] - data[b]) * point;
    }

    #[cfg(feature = "parallel")]
    res.par_iter_mut().enumerate().for_each(|(i, x)| {
        *x = data[i] + (data[i + half_len] - data[i]) * point;
    });

    res
}

macro_rules! zip_self {
    (@ $iter:expr, $step:expr, $skip:expr) => {
        $iter.skip($skip).step_by($step).zip($iter.skip($skip + ($step >> 1)).step_by($step))
    };
    ($iter:expr) => {
        zip_self!(@ $iter, 2, 0)
    };
    ($iter:expr, $step:expr) => {
        zip_self!(@ $iter, $step, 0)
    };
    ($iter:expr, $step:expr, $skip:expr) => {
        zip_self!(@ $iter, $step, $skip)
    };
}

pub(crate) use zip_self;

macro_rules! arc_mpoly {
    ($F:ty, $($elem:expr),* $(,)?) => {
        Arc::new(MultilinearPolynomial::new(vec![$(F::from($elem as u64)),*]))
    };
}

macro_rules! mpoly {
    ($F:ty, $($elem:expr),* $(,)?) => {
        MultilinearPolynomial::<F>::new(vec![$(F::from($elem as u64)),*])
    };
}

macro_rules! arc_mpoly_fr {
    ($($elem:expr),* $(,)?) => {
        Arc::new(MultilinearPolynomial::<F>::new(vec![$($elem),*]))
    };
}

macro_rules! mpoly_fr {
    ($($elem:expr),* $(,)?) => {
        MultilinearPolynomial::<F>::new(vec![$($elem),*])
    };
}

pub(crate) use arc_mpoly;
pub(crate) use arc_mpoly_fr;
pub(crate) use mpoly;
pub(crate) use mpoly_fr;

#[cfg(test)]
mod test {
    use crate::backend::{
        poly::multilinear::{rotation_eval, zip_self, MultilinearPolynomial},
        util::{
            arithmetic::{BooleanHypercube, Field},
            expression::Rotation,
            test::rand_vec,
            Itertools,
        },
    };
    use halo2_curves::bn256::Fr;
    use rand::{rngs::OsRng, RngCore};
    use std::{borrow::Cow, iter};

    #[test]
    fn test_eq_xy() {
        let point = [1, 2, 3, 4].iter().map(|x| Fr::from(*x)).collect_vec();
        let want = (0..16)
            .into_iter()
            .map(|x| {
                let mut res = Fr::one();
                for i in 0..4 {
                    let y = Fr::from((x >> i) & 1);
                    res = res * ((Fr::one() - point[i]) * (Fr::one() - y) + point[i] * y);
                }
                res
            })
            .collect_vec();
        let poly = MultilinearPolynomial::eq_xy(point.as_slice());
        assert_eq!(poly.evals(), want.as_slice());
    }

    #[test]
    fn test_eq_xy_shifted_ab() {
        let point = [1, 2, 3, 4].iter().map(|x| Fr::from(*x)).collect_vec();
        let point_shifted = point[0..point.len() - 1].to_vec();
        for a in 0..=1 as u64 {
            for b in 0..=1 as u64 {
                let f_a = Fr::from(a);
                let f_b = Fr::from(b);
                let point2 = [vec![f_b], point_shifted.clone()].concat();
                let point_l = point[point.len() - 1];
                let want = (0..16)
                    .into_iter()
                    .map(|x| {
                        let mut res = Fr::one();
                        for i in 0..4 {
                            let y = Fr::from((x >> i) & 1);
                            res = res * ((Fr::one() - point2[i]) * (Fr::one() - y) + point2[i] * y);
                        }
                        res * (point_l * f_a + (Fr::one() - point_l) * (Fr::one() - f_a))
                    })
                    .collect_vec();
                let poly =
                    MultilinearPolynomial::eq_xy_shifted_ab(point.as_slice(), a == 1, b == 1);
                assert_eq!(poly.evals(), want.as_slice());
            }
        }
    }

    fn fix_variables<F: Field>(evals: &[F], x: &[F]) -> Vec<F> {
        x.iter()
            .fold(Cow::Borrowed(evals), |evals, x_i| {
                zip_self!(evals.iter())
                    .map(|(eval_0, eval_1)| (*eval_1 - eval_0) * x_i + eval_0)
                    .collect_vec()
                    .into()
            })
            .into_owned()
    }

    #[test]
    fn fix_variable() {
        let rand_x_i = || match OsRng.next_u32() % 3 {
            0 => Fr::zero(),
            1 => Fr::one(),
            2 => Fr::random(OsRng),
            _ => unreachable!(),
        };
        for num_vars in 0..16 {
            let poly = MultilinearPolynomial::rand(num_vars, OsRng);
            for x in (1..=num_vars).map(|n| iter::repeat_with(rand_x_i).take(n).collect_vec()) {
                let mut buf = MultilinearPolynomial::new(vec![Fr::zero(); 1 << (num_vars - 1)]);
                assert_eq!(
                    x.iter()
                        .fold(poly.clone(), |mut poly, x_i| {
                            poly.fix_variable(x_i, &mut buf);
                            poly
                        })
                        .evals(),
                    fix_variables(poly.evals(), &x)
                );
                if x.len() == num_vars {
                    assert_eq!(poly.evaluate(&x), fix_variables(poly.evals(), &x)[0]);
                }
            }
        }
    }

    #[test]
    fn test_fix_variables() {
        let poly = MultilinearPolynomial::new(
            [
                Fr::from(1u64),
                Fr::from(2u64),
                Fr::from(3u64),
                Fr::from(4u64),
            ]
            .to_vec(),
        );
        let x = [Fr::from(2u64), Fr::from(3u64)];
        let want = vec![vec![Fr::from(3u64), Fr::from(5u64)], vec![Fr::from(9u64)]];
        for i in 0..2 {
            let got = poly.fix_variables(&x[..(i + 1)]).into_evals();
            assert_eq!(got, want[i]);
        }
    }

    #[test]
    fn evaluate_for_rotation() {
        let mut rng = OsRng;
        for num_vars in 0..16 {
            let bh = BooleanHypercube::new(num_vars);
            let rotate = |f: &Vec<Fr>| {
                (0..1 << num_vars)
                    .map(|idx| f[bh.rotate(idx, Rotation::next())])
                    .collect_vec()
            };
            let f = rand_vec(1 << num_vars, &mut rng);
            let fs = iter::successors(Some(f), |f| Some(rotate(f)))
                .map(MultilinearPolynomial::new)
                .take(num_vars)
                .collect_vec();
            let x = rand_vec::<Fr>(num_vars, &mut rng);

            for rotation in -(num_vars as i32) + 1..num_vars as i32 {
                let rotation = Rotation(rotation);
                let (f, f_rotated) = if rotation < Rotation::cur() {
                    (fs.last().unwrap(), &fs[fs.len() - rotation.distance() - 1])
                } else {
                    (fs.first().unwrap(), &fs[rotation.distance()])
                };
                assert_eq!(
                    rotation_eval(&x, rotation, &f.evaluate_for_rotation(&x, rotation)),
                    f_rotated.evaluate(&x),
                );
            }
        }
    }
}
