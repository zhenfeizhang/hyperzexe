use crate::fields::{fp, fp12, fp2, FieldExtPoint};
use halo2_base::halo2_proofs::{curves::CurveAffine, group::prime::PrimeCurveAffine};
use halo2_ecc::bigint::CRTInteger;

pub mod final_exp;

type ScalarChip<BF, SF> = fp::FpConfig<BF, SF>;
type ScalarPoint<BF> = CRTInteger<BF>;
type Scalar2Chip<'a, BF, SF, SF2> = fp2::Fp2Chip<'a, BF, ScalarChip<BF, SF>, SF2>;
type Scalar12Chip<'a, BF, SF, SF12> = fp12::Fp12Chip<'a, BF, ScalarChip<BF, SF>, SF12, 9>;

#[cfg(test)]
pub(crate) mod tests;
