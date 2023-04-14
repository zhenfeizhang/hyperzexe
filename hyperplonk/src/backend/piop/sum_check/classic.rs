use crate::{
    backend::{
        piop::sum_check::{SumCheck, VirtualPolynomial},
        poly::multilinear::MultilinearPolynomial,
        util::{
            arithmetic::{BooleanHypercube, Field, PrimeField},
            expression::{Expression, Rotation},
            parallel::par_map_collect,
            start_timer,
            transcript::{FieldTranscriptRead, FieldTranscriptWrite},
            Itertools,
        },
    },
    Error,
};
use num_integer::Integer;
use std::{collections::HashMap, fmt::Debug, marker::PhantomData};

mod coeff;
mod eval;

pub use coeff::CoefficientsProver;
pub use eval::EvaluationsProver;

#[derive(Debug)]
pub struct ProverState<'a, F: Field> {
    num_vars: usize,
    expression: &'a Expression<F>,
    degree: usize,
    sum: F,
    lagranges: HashMap<i32, (usize, F)>,
    identities: Vec<F>,
    eq_xys: Vec<MultilinearPolynomial<F>>,
    eq_xy_shifted_00s: Vec<MultilinearPolynomial<F>>,
    eq_xy_shifted_01s: Vec<MultilinearPolynomial<F>>,
    eq_xy_shifted_10s: Vec<MultilinearPolynomial<F>>,
    eq_xy_shifted_11s: Vec<MultilinearPolynomial<F>>,
    polys: Vec<Vec<MultilinearPolynomial<F>>>,
    challenges: &'a [F],
    buf: MultilinearPolynomial<F>,
    round: usize,
    bh: BooleanHypercube,
}

impl<'a, F: PrimeField> ProverState<'a, F> {
    fn new(num_vars: usize, sum: F, virtual_poly: VirtualPolynomial<'a, F>) -> Self {
        assert!(num_vars > 0 && virtual_poly.expression.max_used_rotation_distance() <= num_vars);
        let bh = BooleanHypercube::new(num_vars);
        let lagranges = {
            // let bh = bh.iter().collect_vec();
            virtual_poly
                .expression
                .used_langrange()
                .into_iter()
                .map(|i| {
                    let b = i.rem_euclid(1 << num_vars as i32) as usize;
                    (i, (b, F::one()))
                })
                .collect()
        };
        let identities = (0..)
            .map(|idx| F::from(idx << num_vars))
            .take(
                virtual_poly
                    .expression
                    .used_identity()
                    .into_iter()
                    .max()
                    .unwrap_or_default()
                    + 1,
            )
            .collect_vec();
        let eq_xys = virtual_poly
            .ys
            .iter()
            .map(|y| MultilinearPolynomial::eq_xy(y))
            .collect_vec();
        let eq_xy_shifted_00s = virtual_poly
            .ys
            .iter()
            .map(|y| MultilinearPolynomial::eq_xy_shifted_ab(y, false, false))
            .collect_vec();
        let eq_xy_shifted_01s = virtual_poly
            .ys
            .iter()
            .map(|y| MultilinearPolynomial::eq_xy_shifted_ab(y, false, true))
            .collect_vec();
        let eq_xy_shifted_10s = virtual_poly
            .ys
            .iter()
            .map(|y| MultilinearPolynomial::eq_xy_shifted_ab(y, true, false))
            .collect_vec();
        let eq_xy_shifted_11s = virtual_poly
            .ys
            .iter()
            .map(|y| MultilinearPolynomial::eq_xy_shifted_ab(y, true, true))
            .collect_vec();
        let polys = virtual_poly
            .polys
            .into_iter()
            .map(|poly| {
                let mut polys = vec![MultilinearPolynomial::zero(); 2 * num_vars];
                polys[num_vars] = MultilinearPolynomial::new(poly.evals().to_vec());
                polys
            })
            .collect_vec();
        Self {
            num_vars,
            expression: virtual_poly.expression,
            degree: virtual_poly.expression.degree(),
            sum,
            lagranges,
            identities,
            eq_xys,
            eq_xy_shifted_00s,
            eq_xy_shifted_01s,
            eq_xy_shifted_10s,
            eq_xy_shifted_11s,
            polys,
            challenges: virtual_poly.challenges,
            buf: MultilinearPolynomial::new(vec![F::zero(); 1 << (num_vars - 1)]),
            round: 0,
            bh,
        }
    }

    fn size(&self) -> usize {
        1 << (self.num_vars - self.round - 1)
    }

    fn next_round(&mut self, sum: F, challenge: &F) {
        self.sum = sum;
        self.lagranges.values_mut().for_each(|(b, value)| {
            if b.is_even() {
                *value *= &(F::one() - challenge);
            } else {
                *value *= challenge;
            }
            *b >>= 1;
        });
        self.identities
            .iter_mut()
            .for_each(|constant| *constant += F::from(1 << self.round) * challenge);
        self.eq_xys
            .iter_mut()
            .for_each(|eq_xy| eq_xy.fix_variable(challenge, &mut self.buf));
        self.eq_xy_shifted_00s
            .iter_mut()
            .for_each(|eq_xy_shift| eq_xy_shift.fix_variable(challenge, &mut self.buf));
        self.eq_xy_shifted_01s
            .iter_mut()
            .for_each(|eq_xy_shift| eq_xy_shift.fix_variable(challenge, &mut self.buf));
        self.eq_xy_shifted_10s
            .iter_mut()
            .for_each(|eq_xy_shift| eq_xy_shift.fix_variable(challenge, &mut self.buf));
        self.eq_xy_shifted_11s
            .iter_mut()
            .for_each(|eq_xy_shift| eq_xy_shift.fix_variable(challenge, &mut self.buf));
        if self.round == 0 {
            let rotation_maps = self
                .expression
                .used_rotation()
                .into_iter()
                .filter_map(|rotation| {
                    (rotation != Rotation::cur())
                        .then(|| (rotation, self.bh.rotation_map(rotation)))
                })
                .collect::<HashMap<_, _>>();
            for query in self.expression.used_query() {
                if query.rotation() != Rotation::cur() {
                    let poly = &self.polys[query.poly()][self.num_vars];
                    let mut rotated = MultilinearPolynomial::new(par_map_collect(
                        &rotation_maps[&query.rotation()],
                        |b| poly[*b],
                    ));
                    rotated.fix_variable(challenge, &mut self.buf);
                    self.polys[query.poly()]
                        [(query.rotation().0 + self.num_vars as i32) as usize] = rotated;
                }
            }
            let size = self.size();
            self.polys.iter_mut().for_each(|polys| {
                let mut output = MultilinearPolynomial::new(vec![F::zero(); size]);
                polys[self.num_vars].fix_variable_into(&mut output, challenge);
                polys[self.num_vars] = output;
            });
        } else {
            self.polys.iter_mut().for_each(|polys| {
                polys.iter_mut().for_each(|poly| {
                    if !poly.is_zero() {
                        poly.fix_variable(challenge, &mut self.buf);
                    }
                });
            });
        }
        self.round += 1;
        self.bh = BooleanHypercube::new(self.num_vars - self.round);
    }

    fn into_evals(self) -> Vec<F> {
        debug_assert_eq!(self.round, self.num_vars);
        self.polys
            .iter()
            .map(|polys| polys[self.num_vars].evals()[0])
            .collect()
    }
}

pub trait ClassicSumCheckProver<F: Field>: Clone + Debug {
    type RoundMessage: ClassicSumCheckRoundMessage<F>;

    fn new(state: &ProverState<F>) -> Self;

    fn prove_round(&self, state: &ProverState<F>) -> Self::RoundMessage;
}

pub trait ClassicSumCheckRoundMessage<F: Field>: Sized + Debug + Clone {
    type Auxiliary: Default;

    fn new(msgs: &[F]) -> Self;

    fn write(&self, transcript: &mut impl FieldTranscriptWrite<F>) -> Result<Vec<F>, Error>;

    fn read(degree: usize, transcript: &mut impl FieldTranscriptRead<F>) -> Result<Self, Error>;

    fn sum(&self) -> F;

    fn auxiliary(_degree: usize) -> Self::Auxiliary {
        Default::default()
    }

    fn evaluate(&self, aux: &Self::Auxiliary, challenge: &F) -> F;

    fn verify_consistency(
        degree: usize,
        mut sum: F,
        msgs: &[Self],
        challenges: &[F],
    ) -> Result<F, Error> {
        let aux = Self::auxiliary(degree);
        for (round, (msg, challenge)) in msgs.iter().zip(challenges.iter()).enumerate() {
            if sum != msg.sum() {
                let msg = if round == 0 {
                    format!("Expect sum {sum:?} but get {:?}", msg.sum())
                } else {
                    format!("Consistency failure at round {round}")
                };
                return Err(Error::InvalidSumcheck(msg));
            }
            sum = msg.evaluate(&aux, challenge);
        }
        Ok(sum)
    }
}

#[derive(Clone, Debug)]
pub struct ClassicSumCheckProof<M>(pub Vec<M>);

#[derive(Clone, Debug)]
pub struct ClassicSumCheck<P>(PhantomData<P>);

impl<F, P> SumCheck<F> for ClassicSumCheck<P>
where
    F: PrimeField,
    P: ClassicSumCheckProver<F>,
{
    type ProverParam = ();
    type VerifierParam = ();
    type Proof = ClassicSumCheckProof<P::RoundMessage>;

    fn prove(
        _: &Self::ProverParam,
        num_vars: usize,
        virtual_poly: VirtualPolynomial<F>,
        sum: F,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(Vec<Vec<F>>, Vec<F>, Vec<F>), Error> {
        let _timer = start_timer(|| {
            let degree = virtual_poly.expression.degree();
            format!("sum_check_prove-{num_vars}-{degree}")
        });

        let mut state = ProverState::new(num_vars, sum, virtual_poly);
        let mut challenges = Vec::with_capacity(num_vars);
        let prover = P::new(&state);
        let aux = P::RoundMessage::auxiliary(state.degree);
        let mut msgs = Vec::with_capacity(num_vars);

        for _ in 0..num_vars {
            let msg = prover.prove_round(&state);
            msgs.push(msg.write(transcript)?);

            let challenge = transcript.squeeze_challenge();
            challenges.push(challenge);

            state.next_round(msg.evaluate(&aux, &challenge), &challenge);
        }

        Ok((msgs, challenges, state.into_evals()))
    }

    fn verify(
        _: &Self::VerifierParam,
        msgs: &[&[F]],
        num_vars: usize,
        degree: usize,
        sum: F,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(F, Vec<F>), Error> {
        let challenges = {
            let mut challenges = Vec::with_capacity(num_vars);
            for i in 0..num_vars {
                P::RoundMessage::new(msgs[i]).write(transcript)?;
                challenges.push(transcript.squeeze_challenge());
            }
            challenges
        };

        let msgs = msgs
            .into_iter()
            .map(|msg| P::RoundMessage::new(msg))
            .collect::<Vec<_>>();
        Ok((
            P::RoundMessage::verify_consistency(degree, sum, &msgs, &challenges)?,
            challenges,
        ))
    }
}
