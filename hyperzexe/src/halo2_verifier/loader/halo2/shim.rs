use crate::halo2_verifier::{
    halo2_proofs::{
        circuit::{Cell, Value},
        ff::PrimeField,
        plonk::Error,
    },
    util::arithmetic::CurveAffine,
};
use std::{fmt::Debug, ops::Deref};

pub trait Context: Debug {
    fn constrain_equal(&mut self, lhs: Cell, rhs: Cell) -> Result<(), Error>;

    fn offset(&self) -> usize;
}

pub trait IntegerInstructions<'a, F: PrimeField>: Clone + Debug {
    type Context: Context;
    type AssignedCell: Clone + Debug;
    type AssignedInteger: Clone + Debug;

    fn assign_scalar(
        &self,
        ctx: &mut Self::Context,
        integer: Value<F>,
    ) -> Result<Self::AssignedInteger, Error>;

    fn assign_constant(
        &self,
        ctx: &mut Self::Context,
        integer: F,
    ) -> Result<Self::AssignedInteger, Error>;

    fn sum_with_coeff_and_const(
        &self,
        ctx: &mut Self::Context,
        values: &[(F, impl Deref<Target = Self::AssignedInteger>)],
        constant: F,
    ) -> Result<Self::AssignedInteger, Error>;

    fn sum_products_with_coeff_and_const(
        &self,
        ctx: &mut Self::Context,
        values: &[(
            F,
            impl Deref<Target = Self::AssignedInteger>,
            impl Deref<Target = Self::AssignedInteger>,
        )],
        constant: F,
    ) -> Result<Self::AssignedInteger, Error>;

    fn sub(
        &self,
        ctx: &mut Self::Context,
        lhs: &Self::AssignedInteger,
        rhs: &Self::AssignedInteger,
    ) -> Result<Self::AssignedInteger, Error>;

    fn neg(
        &self,
        ctx: &mut Self::Context,
        value: &Self::AssignedInteger,
    ) -> Result<Self::AssignedInteger, Error>;

    fn invert(
        &self,
        ctx: &mut Self::Context,
        value: &Self::AssignedInteger,
    ) -> Result<Self::AssignedInteger, Error>;

    fn assert_equal(
        &self,
        ctx: &mut Self::Context,
        lhs: &Self::AssignedInteger,
        rhs: &Self::AssignedInteger,
    ) -> Result<(), Error>;
}

pub trait EccInstructions<'a, C: CurveAffine>: Clone + Debug {
    type Context: Context;
    type ScalarChip: IntegerInstructions<
        'a,
        C::Base,
        Context = Self::Context,
        AssignedCell = Self::AssignedCell,
        AssignedInteger = Self::AssignedScalar,
    >;
    type AssignedCell: Clone + Debug;
    type AssignedScalar: Clone + Debug;
    type AssignedEcPoint: Clone + Debug;

    fn scalar_chip(&self) -> &Self::ScalarChip;

    fn assign_constant(
        &self,
        ctx: &mut Self::Context,
        ec_point: C,
    ) -> Result<Self::AssignedEcPoint, Error>;

    fn assign_point(
        &self,
        ctx: &mut Self::Context,
        ec_point: Value<C>,
    ) -> Result<Self::AssignedEcPoint, Error>;

    fn sum_with_const(
        &self,
        ctx: &mut Self::Context,
        values: &[impl Deref<Target = Self::AssignedEcPoint>],
        constant: C,
    ) -> Result<Self::AssignedEcPoint, Error>;

    fn fixed_base_msm(
        &mut self,
        ctx: &mut Self::Context,
        pairs: &[(impl Deref<Target = Self::AssignedScalar>, C)],
    ) -> Result<Self::AssignedEcPoint, Error>;

    fn variable_base_msm(
        &mut self,
        ctx: &mut Self::Context,
        pairs: &[(
            impl Deref<Target = Self::AssignedScalar>,
            impl Deref<Target = Self::AssignedEcPoint>,
        )],
    ) -> Result<Self::AssignedEcPoint, Error>;

    fn assert_equal(
        &self,
        ctx: &mut Self::Context,
        lhs: &Self::AssignedEcPoint,
        rhs: &Self::AssignedEcPoint,
    ) -> Result<(), Error>;
}

mod halo2_lib {
    use crate::{
        halo2_verifier::{
            halo2_proofs::{
                circuit::{Cell, Value},
                halo2curves::CurveAffineExt,
                plonk::Error,
            },
            loader::halo2::{Context, EccInstructions, IntegerInstructions},
            util::arithmetic::CurveAffine,
        },
        EcPoint, EccChip, FixedEcPoint, ScalarChip, ScalarPoint,
    };
    use halo2_base::{
        self,
        gates::{flex_gate::FlexGateConfig, GateInstructions, RangeInstructions},
        utils::PrimeField,
        AssignedValue,
        QuantumCell::{Constant, Existing, Witness},
    };
    use halo2_proofs::group::prime::PrimeCurveAffine;
    use std::ops::Deref;

    type AssignedScalar<C> = ScalarPoint<C>;
    type AssignedEcPoint<C> = EcPoint<C>;

    impl<'a, F: PrimeField> Context for halo2_base::Context<'a, F> {
        fn constrain_equal(&mut self, lhs: Cell, rhs: Cell) -> Result<(), Error> {
            #[cfg(feature = "halo2-axiom")]
            self.region.constrain_equal(&lhs, &rhs);
            #[cfg(feature = "halo2-pse")]
            self.region.constrain_equal(lhs, rhs)?;
            Ok(())
        }

        fn offset(&self) -> usize {
            unreachable!()
        }
    }

    impl<'a, C: CurveAffine + PrimeCurveAffine> IntegerInstructions<'a, C::Base> for ScalarChip<C> {
        type Context = halo2_base::Context<'a, C::Base>;
        type AssignedCell = AssignedValue<C::Base>;
        type AssignedInteger = AssignedScalar<C>;

        fn assign_scalar(
            &self,
            ctx: &mut Self::Context,
            integer: Value<C::Scalar>,
        ) -> Result<Self::AssignedInteger, Error> {
            Ok(self.assign_region_last(ctx, vec![Witness(integer)], vec![]))
        }

        fn assign_constant(
            &self,
            ctx: &mut Self::Context,
            integer: C::Scalar,
        ) -> Result<Self::AssignedInteger, Error> {
            unimplemented!();
        }

        fn sum_with_coeff_and_const(
            &self,
            ctx: &mut Self::Context,
            values: &[(C::Scalar, impl Deref<Target = Self::AssignedInteger>)],
            constant: C::Scalar,
        ) -> Result<Self::AssignedInteger, Error> {
            unimplemented!();
        }

        fn sum_products_with_coeff_and_const(
            &self,
            ctx: &mut Self::Context,
            values: &[(
                C::Scalar,
                impl Deref<Target = Self::AssignedInteger>,
                impl Deref<Target = Self::AssignedInteger>,
            )],
            constant: C::Scalar,
        ) -> Result<Self::AssignedInteger, Error> {
            unimplemented!();
        }

        fn sub(
            &self,
            ctx: &mut Self::Context,
            a: &Self::AssignedInteger,
            b: &Self::AssignedInteger,
        ) -> Result<Self::AssignedInteger, Error> {
            unimplemented!();
        }

        fn neg(
            &self,
            ctx: &mut Self::Context,
            a: &Self::AssignedInteger,
        ) -> Result<Self::AssignedInteger, Error> {
            unimplemented!();
        }

        fn invert(
            &self,
            ctx: &mut Self::Context,
            a: &Self::AssignedInteger,
        ) -> Result<Self::AssignedInteger, Error> {
            // make sure scalar != 0
            unimplemented!();
        }

        fn assert_equal(
            &self,
            ctx: &mut Self::Context,
            a: &Self::AssignedInteger,
            b: &Self::AssignedInteger,
        ) -> Result<(), Error> {
            unimplemented!();
        }
    }

    impl<'a, C: CurveAffineExt> EccInstructions<'a, C> for EccChip<C::Base>
    where
        C::ScalarExt: PrimeField,
        C::Base: PrimeField,
    {
        type Context = halo2_base::Context<'a, C::Scalar>;
        type ScalarChip = FlexGateConfig<C::Scalar>;
        type AssignedCell = AssignedValue<C::Scalar>;
        type AssignedScalar = AssignedScalar<C>;
        type AssignedEcPoint = AssignedEcPoint<C>;

        fn scalar_chip(&self) -> &Self::ScalarChip {
            self.field_chip.range().gate()
        }

        fn assign_constant(
            &self,
            ctx: &mut Self::Context,
            point: C,
        ) -> Result<Self::AssignedEcPoint, Error> {
            let fixed = FixedEcPoint::<C::Scalar, C>::from_curve(point);
            Ok(FixedEcPoint::assign(fixed, &self.gate, ctx))
        }

        fn assign_point(
            &self,
            ctx: &mut Self::Context,
            point: Value<C>,
        ) -> Result<Self::AssignedEcPoint, Error> {
            let assigned = self.assign_point(ctx, point);
            let is_valid = self.is_on_curve_or_infinity::<C>(ctx, &assigned);
            // Self::ScalarChip.assert_is_const(ctx, &is_valid, C::Scalar::ONE);
            Ok(assigned)
        }

        fn sum_with_const(
            &self,
            ctx: &mut Self::Context,
            values: &[impl Deref<Target = Self::AssignedEcPoint>],
            constant: C,
        ) -> Result<Self::AssignedEcPoint, Error> {
            let constant = if bool::from(constant.is_identity()) {
                None
            } else {
                let constant = EccInstructions::<C>::assign_constant(self, ctx, constant).unwrap();
                Some(constant)
            };
            Ok(self.sum::<C>(ctx, constant.iter().chain(values.iter().map(Deref::deref))))
        }

        fn variable_base_msm(
            &mut self,
            ctx: &mut Self::Context,
            pairs: &[(
                impl Deref<Target = Self::AssignedScalar>,
                impl Deref<Target = Self::AssignedEcPoint>,
            )],
        ) -> Result<Self::AssignedEcPoint, Error> {
            let (scalars, points): (Vec<_>, Vec<_>) = pairs
                .iter()
                .map(|(scalar, point)| (vec![scalar.deref().clone()], point.deref().clone()))
                .unzip();

            Ok(EccChip::<C::Base>::variable_base_msm::<C>(
                self,
                ctx,
                &points,
                &scalars,
                C::Scalar::NUM_BITS as usize,
                4, // empirically clump factor of 4 seems to be best
            ))
        }

        fn fixed_base_msm(
            &mut self,
            ctx: &mut Self::Context,
            pairs: &[(impl Deref<Target = Self::AssignedScalar>, C)],
        ) -> Result<Self::AssignedEcPoint, Error> {
            let (scalars, points): (Vec<_>, Vec<_>) = pairs
                .iter()
                .filter_map(|(scalar, point)| {
                    if point.is_identity().into() {
                        None
                    } else {
                        Some((vec![scalar.deref().clone()], *point))
                    }
                })
                .unzip();

            Ok(EccChip::<C::Base>::fixed_base_msm::<C>(
                self,
                ctx,
                &points,
                &scalars,
                C::Scalar::NUM_BITS as usize,
                0,
                4,
            ))
        }

        fn assert_equal(
            &self,
            ctx: &mut Self::Context,
            a: &Self::AssignedEcPoint,
            b: &Self::AssignedEcPoint,
        ) -> Result<(), Error> {
            self.assert_equal(ctx, a, b);
            Ok(())
        }
    }
}
