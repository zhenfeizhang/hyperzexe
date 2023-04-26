use crate::halo2_verifier::{
    halo2_proofs,
    loader::{
        halo2::{EcPoint, EccInstructions, Halo2Loader, Scalar},
        native::{self, NativeLoader},
        Loader, ScalarLoader,
    },
    util::{
        arithmetic::{fe_to_fe, CurveAffine, PrimeField},
        hash::Poseidon,
        transcript::{Transcript, TranscriptRead, TranscriptWrite},
        Itertools,
    },
    Error,
};
use halo2_proofs::{
    circuit::Value, curves::serde::SerdeObject, ff::FromUniformBytes, transcript::EncodedChallenge,
};
use std::{
    io::{self, Read, Write},
    rc::Rc,
};

/// Encoding that encodes elliptic curve point into native field elements.
pub trait NativeEncoding<'a, C>: EccInstructions<'a, C>
where
    C: CurveAffine,
{
    fn encode(
        &self,
        ctx: &mut Self::Context,
        ec_point: &Self::AssignedEcPoint,
    ) -> Result<Vec<Self::AssignedScalar>, Error>;
}

pub struct PoseidonTranscript<
    C,
    L,
    S,
    const T: usize,
    const RATE: usize,
    const R_F: usize,
    const R_P: usize,
> where
    C: CurveAffine,
    L: Loader<C>,
{
    loader: L,
    stream: S,
    buf: Poseidon<C::Scalar, <L as ScalarLoader<C::Scalar>>::LoadedScalar, T, RATE>,
}

impl<
        'a,
        C,
        R,
        EccChip,
        ScalarChip,
        const T: usize,
        const RATE: usize,
        const R_F: usize,
        const R_P: usize,
    >
    PoseidonTranscript<C, Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>, Value<R>, T, RATE, R_F, R_P>
where
    C: CurveAffine,
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
    R: Read,
    EccChip: NativeEncoding<'a, C>,
{
    pub fn new(loader: &Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>, stream: Value<R>) -> Self {
        let buf = Poseidon::new(loader, R_F, R_P);
        Self {
            loader: loader.clone(),
            stream,
            buf,
        }
    }

    pub fn from_spec(
        loader: &Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>,
        stream: Value<R>,
        spec: crate::halo2_verifier::poseidon::Spec<C::Scalar, T, RATE>,
    ) -> Self {
        let buf = Poseidon::from_spec(loader, spec);
        Self {
            loader: loader.clone(),
            stream,
            buf,
        }
    }

    pub fn new_stream(&mut self, stream: Value<R>) {
        self.buf.clear();
        self.stream = stream;
    }
}

impl<
        'a,
        C,
        R,
        EccChip,
        ScalarChip,
        const T: usize,
        const RATE: usize,
        const R_F: usize,
        const R_P: usize,
    > Transcript<C, Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>>
    for PoseidonTranscript<
        C,
        Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>,
        Value<R>,
        T,
        RATE,
        R_F,
        R_P,
    >
where
    C: CurveAffine,
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
    R: Read,
    EccChip: NativeEncoding<'a, C>,
{
    fn loader(&self) -> &Rc<Halo2Loader<'a, C, EccChip, ScalarChip>> {
        &self.loader
    }

    fn squeeze_challenge(&mut self) -> Scalar<'a, C, EccChip, ScalarChip> {
        self.buf.squeeze()
    }

    fn common_scalar(&mut self, scalar: &Scalar<'a, C, EccChip, ScalarChip>) -> Result<(), Error> {
        self.buf.update(&[scalar.clone()]);
        Ok(())
    }

    fn common_ec_point(
        &mut self,
        ec_point: &EcPoint<'a, C, EccChip, ScalarChip>,
    ) -> Result<(), Error> {
        let encoded = self
            .loader
            .ecc_chip()
            .encode(&mut self.loader.ctx_mut(), &ec_point.assigned())
            .map(|encoded| {
                encoded
                    .into_iter()
                    .map(|encoded| self.loader.scalar_from_assigned(encoded))
                    .collect_vec()
            })
            .map_err(|_| {
                Error::Transcript(
                    io::ErrorKind::Other,
                    "Failed to encode elliptic curve point into native field elements".to_string(),
                )
            })?;
        self.buf.update(&encoded);
        Ok(())
    }
}

impl<
        'a,
        C,
        R,
        EccChip,
        ScalarChip,
        const T: usize,
        const RATE: usize,
        const R_F: usize,
        const R_P: usize,
    > TranscriptRead<C, Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>>
    for PoseidonTranscript<
        C,
        Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>,
        Value<R>,
        T,
        RATE,
        R_F,
        R_P,
    >
where
    C: CurveAffine,
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
    R: Read,
    EccChip: NativeEncoding<'a, C>,
{
    fn read_scalar(&mut self) -> Result<Scalar<'a, C, EccChip, ScalarChip>, Error> {
        let scalar = self.stream.as_mut().and_then(|stream| {
            let mut data = <C::Scalar as PrimeField>::Repr::default();
            if stream.read_exact(data.as_mut()).is_err() {
                return Value::unknown();
            }
            Option::<C::Scalar>::from(C::Scalar::from_repr(data))
                .map(Value::known)
                .unwrap_or_else(Value::unknown)
        });
        let scalar = self.loader.assign_scalar(scalar);
        self.common_scalar(&scalar)?;
        Ok(scalar)
    }

    fn read_ec_point(&mut self) -> Result<EcPoint<'a, C, EccChip, ScalarChip>, Error> {
        let ec_point = self.stream.as_mut().and_then(|stream| {
            let mut compressed = C::Repr::default();
            if stream.read_exact(compressed.as_mut()).is_err() {
                return Value::unknown();
            }
            Option::<C>::from(C::from_bytes(&compressed))
                .map(Value::known)
                .unwrap_or_else(Value::unknown)
        });
        let ec_point = self.loader.assign_ec_point(ec_point);
        self.common_ec_point(&ec_point)?;
        Ok(ec_point)
    }
}

impl<C: CurveAffine, S, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    PoseidonTranscript<C, NativeLoader, S, T, RATE, R_F, R_P>
where
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
{
    pub fn new(stream: S) -> Self {
        Self {
            loader: NativeLoader,
            stream,
            buf: Poseidon::new(&NativeLoader, R_F, R_P),
        }
    }

    pub fn from_spec(
        stream: S,
        spec: crate::halo2_verifier::poseidon::Spec<C::Scalar, T, RATE>,
    ) -> Self {
        Self {
            loader: NativeLoader,
            stream,
            buf: Poseidon::from_spec(&NativeLoader, spec),
        }
    }

    pub fn new_stream(&mut self, stream: S) {
        self.buf.clear();
        self.stream = stream;
    }
}

impl<C: CurveAffine, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    PoseidonTranscript<C, NativeLoader, Vec<u8>, T, RATE, R_F, R_P>
where
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
{
    pub fn clear(&mut self) {
        self.buf.clear();
        self.stream.clear();
    }
}

impl<C: CurveAffine, S, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    Transcript<C, NativeLoader> for PoseidonTranscript<C, NativeLoader, S, T, RATE, R_F, R_P>
where
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
{
    fn loader(&self) -> &NativeLoader {
        &native::LOADER
    }

    fn squeeze_challenge(&mut self) -> C::Scalar {
        self.buf.squeeze()
    }

    fn common_scalar(&mut self, scalar: &C::Scalar) -> Result<(), Error> {
        self.buf.update(&[*scalar]);
        Ok(())
    }

    fn common_ec_point(&mut self, ec_point: &C) -> Result<(), Error> {
        let encoded: Vec<_> = Option::from(ec_point.coordinates().map(|coordinates| {
            [coordinates.x(), coordinates.y()]
                .into_iter()
                .cloned()
                .map(fe_to_fe)
                .collect_vec()
        }))
        .ok_or_else(|| {
            Error::Transcript(
                io::ErrorKind::Other,
                "Invalid elliptic curve point encoding in proof".to_string(),
            )
        })?;
        self.buf.update(&encoded);
        Ok(())
    }
}

impl<C, R, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    TranscriptRead<C, NativeLoader> for PoseidonTranscript<C, NativeLoader, R, T, RATE, R_F, R_P>
where
    C: CurveAffine,
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
    R: Read,
{
    fn read_scalar(&mut self) -> Result<C::Scalar, Error> {
        let mut data = <C::Scalar as PrimeField>::Repr::default();
        self.stream
            .read_exact(data.as_mut())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        let scalar = C::Scalar::from_repr_vartime(data).ok_or_else(|| {
            Error::Transcript(
                io::ErrorKind::Other,
                "Invalid scalar encoding in proof".to_string(),
            )
        })?;
        self.common_scalar(&scalar)?;
        Ok(scalar)
    }

    fn read_ec_point(&mut self) -> Result<C, Error> {
        let mut data = C::Repr::default();
        self.stream
            .read_exact(data.as_mut())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        let ec_point = Option::<C>::from(C::from_bytes(&data)).ok_or_else(|| {
            Error::Transcript(
                io::ErrorKind::Other,
                "Invalid elliptic curve point encoding in proof".to_string(),
            )
        })?;
        self.common_ec_point(&ec_point)?;
        Ok(ec_point)
    }
}

impl<C, W, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    PoseidonTranscript<C, NativeLoader, W, T, RATE, R_F, R_P>
where
    C: CurveAffine,
    W: Write,
{
    pub fn stream_mut(&mut self) -> &mut W {
        &mut self.stream
    }

    pub fn finalize(self) -> W {
        self.stream
    }
}

impl<C, W, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize> TranscriptWrite<C>
    for PoseidonTranscript<C, NativeLoader, W, T, RATE, R_F, R_P>
where
    C: CurveAffine,
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
    W: Write,
{
    fn write_scalar(&mut self, scalar: C::Scalar) -> Result<(), Error> {
        self.common_scalar(&scalar)?;
        let data = scalar.to_repr();
        self.stream_mut().write_all(data.as_ref()).map_err(|err| {
            Error::Transcript(
                err.kind(),
                "Failed to write scalar to transcript".to_string(),
            )
        })
    }

    fn write_ec_point(&mut self, ec_point: C) -> Result<(), Error> {
        self.common_ec_point(&ec_point)?;
        let data = ec_point.to_bytes();
        self.stream_mut().write_all(data.as_ref()).map_err(|err| {
            Error::Transcript(
                err.kind(),
                "Failed to write elliptic curve to transcript".to_string(),
            )
        })
    }
}

pub struct ChallengeScalar<C: CurveAffine>(C::Scalar);

impl<C: CurveAffine> EncodedChallenge<C> for ChallengeScalar<C> {
    type Input = C::Scalar;

    fn new(challenge_input: &C::Scalar) -> Self {
        ChallengeScalar(*challenge_input)
    }

    fn get_scalar(&self) -> C::Scalar {
        self.0
    }
}

impl<C: CurveAffine, S, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    halo2_proofs::transcript::Transcript<C, ChallengeScalar<C>>
    for PoseidonTranscript<C, NativeLoader, S, T, RATE, R_F, R_P>
where
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
{
    fn squeeze_challenge(&mut self) -> ChallengeScalar<C> {
        ChallengeScalar::new(&Transcript::squeeze_challenge(self))
    }

    fn common_point(&mut self, ec_point: C) -> io::Result<()> {
        match Transcript::common_ec_point(self, &ec_point) {
            Err(Error::Transcript(kind, msg)) => Err(io::Error::new(kind, msg)),
            Err(_) => unreachable!(),
            _ => Ok(()),
        }
    }

    fn common_scalar(&mut self, scalar: C::Scalar) -> io::Result<()> {
        match Transcript::common_scalar(self, &scalar) {
            Err(Error::Transcript(kind, msg)) => Err(io::Error::new(kind, msg)),
            Err(_) => unreachable!(),
            _ => Ok(()),
        }
    }
}

impl<C, R, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    halo2_proofs::transcript::TranscriptRead<C, ChallengeScalar<C>>
    for PoseidonTranscript<C, NativeLoader, R, T, RATE, R_F, R_P>
where
    C: CurveAffine,
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
    R: Read,
{
    fn read_point(&mut self) -> io::Result<C> {
        match TranscriptRead::read_ec_point(self) {
            Err(Error::Transcript(kind, msg)) => Err(io::Error::new(kind, msg)),
            Err(_) => unreachable!(),
            Ok(value) => Ok(value),
        }
    }

    fn read_scalar(&mut self) -> io::Result<C::Scalar> {
        match TranscriptRead::read_scalar(self) {
            Err(Error::Transcript(kind, msg)) => Err(io::Error::new(kind, msg)),
            Err(_) => unreachable!(),
            Ok(value) => Ok(value),
        }
    }
}

impl<C, R, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    halo2_proofs::transcript::TranscriptReadBuffer<R, C, ChallengeScalar<C>>
    for PoseidonTranscript<C, NativeLoader, R, T, RATE, R_F, R_P>
where
    C: CurveAffine,
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
    R: Read,
{
    fn init(reader: R) -> Self {
        Self::new(reader)
    }
}

impl<C, W, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    halo2_proofs::transcript::TranscriptWrite<C, ChallengeScalar<C>>
    for PoseidonTranscript<C, NativeLoader, W, T, RATE, R_F, R_P>
where
    C: CurveAffine,
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
    W: Write,
{
    fn write_point(&mut self, ec_point: C) -> io::Result<()> {
        halo2_proofs::transcript::Transcript::<C, ChallengeScalar<C>>::common_point(
            self, ec_point,
        )?;
        let data = ec_point.to_bytes();
        self.stream_mut().write_all(data.as_ref())
    }

    fn write_scalar(&mut self, scalar: C::Scalar) -> io::Result<()> {
        halo2_proofs::transcript::Transcript::<C, ChallengeScalar<C>>::common_scalar(self, scalar)?;
        let data = scalar.to_repr();
        self.stream_mut().write_all(data.as_ref())
    }
}

impl<C, W, const T: usize, const RATE: usize, const R_F: usize, const R_P: usize>
    halo2_proofs::transcript::TranscriptWriterBuffer<W, C, ChallengeScalar<C>>
    for PoseidonTranscript<C, NativeLoader, W, T, RATE, R_F, R_P>
where
    C: CurveAffine,
    C::ScalarExt: SerdeObject + FromUniformBytes<64>,
    W: Write,
{
    fn init(writer: W) -> Self {
        Self::new(writer)
    }

    fn finalize(self) -> W {
        self.finalize()
    }
}

mod halo2_lib {
    use crate::{
        halo2_verifier::{
            halo2_curves::CurveAffineExt, system::halo2::transcript::halo2::NativeEncoding,
        },
        EccChip,
    };
    use halo2_base::utils::PrimeField;

    impl<'a, C: CurveAffineExt> NativeEncoding<'a, C> for EccChip<C::Base>
    where
        C::Scalar: PrimeField,
        C::Base: PrimeField,
    {
        fn encode(
            &self,
            _: &mut Self::Context,
            ec_point: &Self::AssignedEcPoint,
        ) -> Result<Vec<Self::AssignedScalar>, crate::halo2_verifier::Error> {
            Ok(vec![
                ec_point.x().native().clone(),
                ec_point.y().native().clone(),
            ])
        }
    }
}
