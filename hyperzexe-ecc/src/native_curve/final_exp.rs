use super::{Scalar12Chip, Scalar2Chip, ScalarChip, ScalarPoint};
use crate::{
    ecc::get_naf,
    fields::{fp12::mul_no_carry_w6, FieldChip, FieldExtPoint},
    halo2_proofs::arithmetic::Field,
};
use halo2_base::{
    gates::GateInstructions,
    halo2_proofs::{
        curves::{CurveAffine, CurveAffineExt},
        ff::PrimeField,
        group::prime::{PrimeCurve, PrimeCurveAffine},
    },
    utils::{fe_to_biguint, modulus, PrimeField},
    Context,
    QuantumCell::{Constant, Existing},
};
use num_bigint::BigUint;

const XI_0: i64 = 9;

impl<'a, BF, SF, SF12> Scalar12Chip<'a, BF, SF, SF12>
where
    BF: PrimeField,
    SF: PrimeField,
    SF12: Field,
{
    // computes a ** (p ** power)
    // only works for p = 3 (mod 4) and p = 1 (mod 6)
    pub fn frobenius_map(
        &self,
        ctx: &mut Context<'_, C::Base>,
        a: &<Self as FieldChip<C::Base>>::FieldPoint,
        power: usize,
    ) -> <Self as FieldChip<C::Base>>::FieldPoint {
        unimplemented!();
    }

    // exp is in little-endian
    pub fn pow(
        &self,
        ctx: &mut Context<'_, C::Base>,
        a: &<Self as FieldChip<C::Base>>::FieldPoint,
        exp: Vec<u64>,
    ) -> <Self as FieldChip<C::Base>>::FieldPoint {
        unimplemented!();
    }

    // assume input is an element of Fp12 in the cyclotomic subgroup GΦ₁₂
    // A cyclotomic group is a subgroup of Fp^n defined by
    //   GΦₙ(p) = {α ∈ Fpⁿ : α^{Φₙ(p)} = 1}

    // below we implement compression and decompression for an element  GΦ₁₂ following Theorem 3.1 of https://eprint.iacr.org/2010/542.pdf
    // Fp4 = Fp2(w^3) where (w^3)^2 = XI_0 +u
    // Fp12 = Fp4(w) where w^3 = w^3

    /// in = g0 + g2 w + g4 w^2 + g1 w^3 + g3 w^4 + g5 w^5 where g_i = g_i0 +
    /// g_i1 * u are elements of Fp2 out = Compress(in) = [ g2, g3, g4, g5 ]
    pub fn cyclotomic_compress(
        &self,
        a: &FieldExtPoint<ScalarPoint<C>>,
    ) -> Vec<FieldExtPoint<ScalarPoint<C>>> {
        let g2 = FieldExtPoint::construct(vec![a.coeffs[1].clone(), a.coeffs[1 + 6].clone()]);
        let g3 = FieldExtPoint::construct(vec![a.coeffs[4].clone(), a.coeffs[4 + 6].clone()]);
        let g4 = FieldExtPoint::construct(vec![a.coeffs[2].clone(), a.coeffs[2 + 6].clone()]);
        let g5 = FieldExtPoint::construct(vec![a.coeffs[5].clone(), a.coeffs[5 + 6].clone()]);
        vec![g2, g3, g4, g5]
    }

    /// Input:
    /// * `compression = [g2, g3, g4, g5]` where g_i are proper elements of Fp2
    /// Output:
    /// * `Decompress(compression) = g0 + g2 w + g4 w^2 + g1 w^3 + g3 w^4 + g5
    ///   w^5` where
    /// * All elements of output are proper elements of Fp2 and: c = XI0 + u if
    ///   g2 != 0: g1 = (g5^2 * c + 3 g4^2 - 2 g3)/(4g2) g0 = (2 g1^2 + g2 * g5
    ///   - 3 g3*g4) * c + 1 if g2 = 0: g1 = (2 g4 * g5)/g3 g0 = (2 g1^2 - 3 g3
    ///   * g4) * c + 1
    pub fn cyclotomic_decompress(
        &self,
        ctx: &mut Context<'_, C::Base>,
        compression: Vec<FieldExtPoint<ScalarPoint<C>>>,
    ) -> FieldExtPoint<ScalarPoint<C>> {
        unimplemented!();
    }

    // input is [g2, g3, g4, g5] = C(g) in compressed format of
    // `cyclotomic_compress` assume all inputs are proper Fp2 elements
    // output is C(g^2) = [h2, h3, h4, h5] computed using Theorem 3.2 of https://eprint.iacr.org/2010/542.pdf
    // all output elements are proper Fp2 elements (with carry)
    //  c = XI_0 + u
    //  h2 = 2(g2 + 3*c*B_45)
    //  h3 = 3(A_45 - (c+1)B_45) - 2g3
    //  h4 = 3(A_23 - (c+1)B_23) - 2g4
    //  h5 = 2(g5 + 3B_23)
    //  A_ij = (g_i + g_j)(g_i + c g_j)
    //  B_ij = g_i g_j

    pub fn cyclotomic_square(
        &self,
        ctx: &mut Context<'_, C::Base>,
        compression: &[FieldExtPoint<ScalarPoint<C>>],
    ) -> Vec<FieldExtPoint<ScalarPoint<C>>> {
        unimplemented!();
    }

    // exp is in little-endian
    pub fn cyclotomic_pow(
        &self,
        ctx: &mut Context<'_, C::Base>,
        a: FieldExtPoint<ScalarPoint<C>>,
        exp: Vec<u64>,
    ) -> FieldExtPoint<ScalarPoint<C>> {
        unimplemented!();
    }

    #[allow(non_snake_case)]
    // use equation for (p^4 - p^2 + 1)/r in Section 5 of https://eprint.iacr.org/2008/490.pdf for BN curves
    pub fn hard_part_BN(
        &self,
        ctx: &mut Context<'_, C::Base>,
        m: <Self as FieldChip<C::Base>>::FieldPoint,
    ) -> <Self as FieldChip<C::Base>>::FieldPoint {
        unimplemented!();
    }

    // out = in^{ (q^6 - 1)*(q^2 + 1) }
    pub fn easy_part(
        &self,
        ctx: &mut Context<'_, C::Base>,
        a: &<Self as FieldChip<C::Base>>::FieldPoint,
    ) -> <Self as FieldChip<C::Base>>::FieldPoint {
        unimplemented!();
    }

    // out = in^{(q^12 - 1)/r}
    pub fn final_exp(
        &self,
        ctx: &mut Context<'_, C::Base>,
        a: &<Self as FieldChip<C::Base>>::FieldPoint,
    ) -> <Self as FieldChip<C::Base>>::FieldPoint {
        let f0 = self.easy_part(ctx, a);
        let f = self.hard_part_BN(ctx, f0);
        f
    }
}
