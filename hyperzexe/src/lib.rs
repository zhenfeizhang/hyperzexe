//! Main module for the Hyperzexe

use halo2_base::gates::flex_gate::FlexGateConfig;
use halo2_ecc::fields::FieldChip;
use halo2_proofs::curves::CurveAffine;

mod halo2_verifier;
mod hyperplonk_verifier;

// may use other gates in the future
pub(crate) type EccChip<C: CurveAffine> =
    hyperzexe_ecc::native_ecc::EccChip<C::Base, FlexGateConfig<C::Base>>;
pub(crate) type ScalarChip<C: CurveAffine> = halo2_ecc::fields::fp::FpConfig<C::Base, C::Scalar>;
pub(crate) type EcPoint<C: CurveAffine> = hyperzexe_ecc::native_ecc::EcPoint<C::Base>;
pub(crate) type FixedEcPoint<C: CurveAffine> =
    hyperzexe_ecc::native_ecc::fixed_base::FixedEcPoint<C>;
pub(crate) type ScalarPoint<C: CurveAffine> = <ScalarChip<C> as FieldChip<C::Base>>::FieldPoint;
