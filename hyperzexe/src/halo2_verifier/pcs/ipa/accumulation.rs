use crate::{
    loader::{native::NativeLoader, LoadedScalar, Loader},
    pcs::{
        ipa::{
            h_coeffs, h_eval, Ipa, IpaAccumulator, IpaProof, IpaProvingKey, IpaSuccinctVerifyingKey,
        },
        AccumulationScheme, AccumulationSchemeProver, PolynomialCommitmentScheme,
    },
    util::{
        arithmetic::{Curve, CurveAffine, Field},
        msm::Msm,
        poly::Polynomial,
        transcript::{TranscriptRead, TranscriptWrite},
        Itertools,
    },
    Error,
};
use rand::Rng;
use std::{array, iter, marker::PhantomData};

#[derive(Clone, Debug)]
pub struct IpaAs<PCS>(PhantomData<PCS>);

impl<C, L, PCS> AccumulationScheme<C, L, PCS> for IpaAs<PCS>
where
    C: CurveAffine,
    L: Loader<C>,
    PCS: PolynomialCommitmentScheme<C, L, Accumulator = IpaAccumulator<C, L>>,
{
    type VerifyingKey = IpaSuccinctVerifyingKey<C>;
    type Proof = IpaAsProof<C, L, PCS>;

    fn read_proof<T>(
        vk: &Self::VerifyingKey,
        instances: &[PCS::Accumulator],
        transcript: &mut T,
    ) -> Result<Self::Proof, Error>
    where
        T: TranscriptRead<C, L>,
    {
        IpaAsProof::read(vk, instances, transcript)
    }

    fn verify(
        vk: &Self::VerifyingKey,
        instances: &[PCS::Accumulator],
        proof: &Self::Proof,
    ) -> Result<PCS::Accumulator, Error> {
        let loader = proof.z.loader();
        let s = vk.s.as_ref().map(|s| loader.ec_point_load_const(s));

        let (u, h) = instances
            .iter()
            .map(|IpaAccumulator { u, xi }| (u.clone(), h_eval(xi, &proof.z)))
            .chain(proof.a_b_u.as_ref().map(|(a, b, u)| (u.clone(), a.clone() * &proof.z + b)))
            .unzip::<_, _, Vec<_>, Vec<_>>();
        let powers_of_alpha = proof.alpha.powers(u.len());

        let mut c = powers_of_alpha
            .iter()
            .zip(u.iter())
            .map(|(power_of_alpha, u)| Msm::<C, L>::base(u) * power_of_alpha)
            .sum::<Msm<_, _>>();
        if let Some(omega) = proof.omega.as_ref() {
            c += Msm::base(s.as_ref().unwrap()) * omega;
        }
        let v = loader.sum_products(&powers_of_alpha.iter().zip(h.iter()).collect_vec());

        Ipa::<C, ()>::succinct_verify(vk, &c, &proof.z, &v, &proof.ipa)
    }
}

#[derive(Clone, Debug)]
pub struct IpaAsProof<C, L, PCS>
where
    C: CurveAffine,
    L: Loader<C>,
    PCS: PolynomialCommitmentScheme<C, L, Accumulator = IpaAccumulator<C, L>>,
{
    a_b_u: Option<(L::LoadedScalar, L::LoadedScalar, L::LoadedEcPoint)>,
    omega: Option<L::LoadedScalar>,
    alpha: L::LoadedScalar,
    z: L::LoadedScalar,
    ipa: IpaProof<C, L>,
    _marker: PhantomData<PCS>,
}

impl<C, L, PCS> IpaAsProof<C, L, PCS>
where
    C: CurveAffine,
    L: Loader<C>,
    PCS: PolynomialCommitmentScheme<C, L, Accumulator = IpaAccumulator<C, L>>,
{
    fn read<T>(
        vk: &IpaSuccinctVerifyingKey<C>,
        instances: &[PCS::Accumulator],
        transcript: &mut T,
    ) -> Result<Self, Error>
    where
        T: TranscriptRead<C, L>,
    {
        assert!(instances.len() > 1);

        let a_b_u = vk
            .zk()
            .then(|| {
                let a = transcript.read_scalar()?;
                let b = transcript.read_scalar()?;
                let u = transcript.read_ec_point()?;
                Ok((a, b, u))
            })
            .transpose()?;
        let omega = vk
            .zk()
            .then(|| {
                let omega = transcript.read_scalar()?;
                Ok(omega)
            })
            .transpose()?;

        for accumulator in instances {
            for xi in accumulator.xi.iter() {
                transcript.common_scalar(xi)?;
            }
            transcript.common_ec_point(&accumulator.u)?;
        }

        let alpha = transcript.squeeze_challenge();
        let z = transcript.squeeze_challenge();

        let ipa = IpaProof::read(vk, transcript)?;

        Ok(Self { a_b_u, omega, alpha, z, ipa, _marker: PhantomData })
    }
}

impl<C, PCS> AccumulationSchemeProver<C, PCS> for IpaAs<PCS>
where
    C: CurveAffine,
    PCS: PolynomialCommitmentScheme<C, NativeLoader, Accumulator = IpaAccumulator<C, NativeLoader>>,
{
    type ProvingKey = IpaProvingKey<C>;

    fn create_proof<T, R>(
        pk: &Self::ProvingKey,
        instances: &[PCS::Accumulator],
        transcript: &mut T,
        mut rng: R,
    ) -> Result<PCS::Accumulator, Error>
    where
        T: TranscriptWrite<C>,
        R: Rng,
    {
        assert!(instances.len() > 1);

        let a_b_u = pk
            .zk()
            .then(|| {
                let [a, b] = array::from_fn(|_| C::Scalar::random(&mut rng));
                let u = (pk.g[1] * a + pk.g[0] * b).to_affine();
                transcript.write_scalar(a)?;
                transcript.write_scalar(b)?;
                transcript.write_ec_point(u)?;
                Ok((a, b, u))
            })
            .transpose()?;
        let omega = pk
            .zk()
            .then(|| {
                let omega = C::Scalar::random(&mut rng);
                transcript.write_scalar(omega)?;
                Ok(omega)
            })
            .transpose()?;

        for accumulator in instances {
            for xi in accumulator.xi.iter() {
                transcript.common_scalar(xi)?;
            }
            transcript.common_ec_point(&accumulator.u)?;
        }

        let alpha = transcript.squeeze_challenge();
        let z = transcript.squeeze_challenge();

        let (u, h) = instances
            .iter()
            .map(|IpaAccumulator { u, xi }| (*u, h_coeffs(xi, C::Scalar::ONE)))
            .chain(a_b_u.map(|(a, b, u)| {
                (
                    u,
                    iter::empty()
                        .chain([b, a])
                        .chain(iter::repeat_with(C::Scalar::zero).take(pk.domain.n - 2))
                        .collect(),
                )
            }))
            .unzip::<_, _, Vec<_>, Vec<_>>();
        let powers_of_alpha = alpha.powers(u.len());

        let h = powers_of_alpha
            .into_iter()
            .zip(h.into_iter().map(Polynomial::new))
            .map(|(power_of_alpha, h)| h * power_of_alpha)
            .sum::<Polynomial<_>>();

        Ipa::<C, ()>::create_proof(pk, &h.to_vec(), &z, omega.as_ref(), transcript, &mut rng)
    }
}

#[cfg(test)]
mod test {
    use crate::halo2_curves::pasta::pallas;
    use crate::halo2_proofs::transcript::{
        Blake2bRead, Blake2bWrite, TranscriptReadBuffer, TranscriptWriterBuffer,
    };
    use crate::{
        pcs::{
            ipa::{self, IpaProvingKey},
            AccumulationScheme, AccumulationSchemeProver, Decider,
        },
        util::{arithmetic::Field, msm::Msm, poly::Polynomial, Itertools},
    };
    use rand::rngs::OsRng;
    use std::iter;

    #[test]
    fn test_ipa_as() {
        type Ipa = ipa::Ipa<pallas::Affine, ()>;
        type IpaAs = ipa::IpaAs<Ipa>;

        let k = 10;
        let zk = true;
        let mut rng = OsRng;

        let pk = IpaProvingKey::<pallas::Affine>::rand(k, zk, &mut rng);
        let accumulators = iter::repeat_with(|| {
            let (c, z, v, proof) = {
                let p = Polynomial::<pallas::Scalar>::rand(pk.domain.n, &mut rng);
                let omega = pk.zk().then(|| pallas::Scalar::random(&mut rng));
                let c = pk.commit(&p, omega);
                let z = pallas::Scalar::random(&mut rng);
                let v = p.evaluate(z);
                let mut transcript = Blake2bWrite::init(Vec::new());
                Ipa::create_proof(&pk, &p[..], &z, omega.as_ref(), &mut transcript, &mut rng)
                    .unwrap();
                (c, z, v, transcript.finalize())
            };

            let svk = pk.svk();
            let accumulator = {
                let mut transcript = Blake2bRead::init(proof.as_slice());
                let proof = Ipa::read_proof(&svk, &mut transcript).unwrap();
                Ipa::succinct_verify(&svk, &Msm::base(&c), &z, &v, &proof).unwrap()
            };

            accumulator
        })
        .take(10)
        .collect_vec();

        let proof = {
            let apk = pk.clone();
            let mut transcript = Blake2bWrite::init(Vec::new());
            IpaAs::create_proof(&apk, &accumulators, &mut transcript, &mut rng).unwrap();
            transcript.finalize()
        };

        let accumulator = {
            let avk = pk.svk();
            let mut transcript = Blake2bRead::init(proof.as_slice());
            let proof = IpaAs::read_proof(&avk, &accumulators, &mut transcript).unwrap();
            IpaAs::verify(&avk, &accumulators, &proof).unwrap()
        };

        let dk = pk.dk();
        assert!(Ipa::decide(&dk, accumulator));
    }
}
