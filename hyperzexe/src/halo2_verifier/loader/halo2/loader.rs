use crate::halo2_verifier::{
    halo2_proofs::circuit,
    loader::{
        halo2::shim::{EccInstructions, IntegerInstructions},
        EcPointLoader, LoadedEcPoint, LoadedScalar, Loader, ScalarLoader,
    },
    util::{
        arithmetic::{CurveAffine, Field, FieldOps},
        Itertools,
    },
};
use std::{
    cell::{Ref, RefCell, RefMut},
    fmt::{self, Debug},
    marker::PhantomData,
    ops::{Add, AddAssign, Deref, Mul, MulAssign, Neg, Sub, SubAssign},
    rc::Rc,
};

#[derive(Debug)]
pub struct Halo2Loader<
    'a,
    C: CurveAffine,
    EccChip: EccInstructions<'a, C>,
    ScalarChip: IntegerInstructions<'a, C::Base>,
> {
    ecc_chip: RefCell<EccChip>,
    scalar_chip: RefCell<ScalarChip>,
    ctx: RefCell<EccChip::Context>,
    num_scalar: RefCell<usize>,
    num_ec_point: RefCell<usize>,
    _marker: PhantomData<C>,
    #[cfg(test)]
    row_meterings: RefCell<Vec<(String, usize)>>,
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Halo2Loader<'a, C, EccChip, ScalarChip>
{
    pub fn new(ecc_chip: EccChip, scalar_chip: ScalarChip, ctx: EccChip::Context) -> Rc<Self> {
        Rc::new(Self {
            ecc_chip: RefCell::new(ecc_chip),
            scalar_chip: RefCell::new(scalar_chip),
            ctx: RefCell::new(ctx),
            num_scalar: RefCell::default(),
            num_ec_point: RefCell::default(),
            #[cfg(test)]
            row_meterings: RefCell::default(),
            _marker: PhantomData,
        })
    }

    pub fn into_ctx(self) -> EccChip::Context {
        self.ctx.into_inner()
    }

    pub fn ecc_chip(&self) -> Ref<EccChip> {
        self.ecc_chip.borrow()
    }

    pub fn scalar_chip(&self) -> Ref<EccChip::ScalarChip> {
        self.scalar_chip.borrow()
    }

    pub fn ctx(&self) -> Ref<EccChip::Context> {
        self.ctx.borrow()
    }

    pub fn ctx_mut(&self) -> RefMut<'_, EccChip::Context> {
        self.ctx.borrow_mut()
    }

    fn assign_const_scalar(self: &Rc<Self>, constant: C::Scalar) -> EccChip::AssignedScalar {
        self.scalar_chip()
            .assign_constant(&mut self.ctx_mut(), constant)
            .unwrap()
    }

    pub fn assign_scalar(
        self: &Rc<Self>,
        scalar: circuit::Value<C::Scalar>,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        let assigned = self
            .scalar_chip()
            .assign_scalar(&mut self.ctx_mut(), scalar)
            .unwrap();
        self.scalar_from_assigned(assigned)
    }

    pub fn scalar_from_assigned(
        self: &Rc<Self>,
        assigned: EccChip::AssignedScalar,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        self.scalar(Value::Assigned(assigned))
    }

    fn scalar(
        self: &Rc<Self>,
        value: Value<C::Scalar, EccChip::AssignedScalar>,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        let index = *self.num_scalar.borrow();
        *self.num_scalar.borrow_mut() += 1;
        Scalar {
            loader: self.clone(),
            index,
            value: value.into(),
        }
    }

    fn assign_const_ec_point(self: &Rc<Self>, constant: C) -> EccChip::AssignedEcPoint {
        self.ecc_chip()
            .assign_constant(&mut self.ctx_mut(), constant)
            .unwrap()
    }

    pub fn assign_ec_point(
        self: &Rc<Self>,
        ec_point: circuit::Value<C>,
    ) -> EcPoint<'a, C, EccChip, ScalarChip> {
        let assigned = self
            .ecc_chip()
            .assign_point(&mut self.ctx_mut(), ec_point)
            .unwrap();
        self.ec_point_from_assigned(assigned)
    }

    pub fn ec_point_from_assigned(
        self: &Rc<Self>,
        assigned: EccChip::AssignedEcPoint,
    ) -> EcPoint<'a, C, EccChip, ScalarChip> {
        self.ec_point(Value::Assigned(assigned))
    }

    fn ec_point(
        self: &Rc<Self>,
        value: Value<C, EccChip::AssignedEcPoint>,
    ) -> EcPoint<'a, C, EccChip, ScalarChip> {
        let index = *self.num_ec_point.borrow();
        *self.num_ec_point.borrow_mut() += 1;
        EcPoint {
            loader: self.clone(),
            index,
            value: value.into(),
        }
    }

    fn add(
        self: &Rc<Self>,
        lhs: &Scalar<'a, C, EccChip, ScalarChip>,
        rhs: &Scalar<'a, C, EccChip, ScalarChip>,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        let output = match (lhs.value().deref(), rhs.value().deref()) {
            (Value::Constant(lhs), Value::Constant(rhs)) => Value::Constant(*lhs + rhs),
            (Value::Assigned(assigned), Value::Constant(constant))
            | (Value::Constant(constant), Value::Assigned(assigned)) => self
                .scalar_chip()
                .sum_with_coeff_and_const(
                    &mut self.ctx_mut(),
                    &[(C::Scalar::ONE, assigned)],
                    *constant,
                )
                .map(Value::Assigned)
                .unwrap(),
            (Value::Assigned(lhs), Value::Assigned(rhs)) => self
                .scalar_chip()
                .sum_with_coeff_and_const(
                    &mut self.ctx_mut(),
                    &[(C::Scalar::ONE, lhs), (C::Scalar::ONE, rhs)],
                    C::Scalar::ZERO,
                )
                .map(Value::Assigned)
                .unwrap(),
        };
        self.scalar(output)
    }

    fn sub(
        self: &Rc<Self>,
        lhs: &Scalar<'a, C, EccChip, ScalarChip>,
        rhs: &Scalar<'a, C, EccChip, ScalarChip>,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        let output = match (lhs.value().deref(), rhs.value().deref()) {
            (Value::Constant(lhs), Value::Constant(rhs)) => Value::Constant(*lhs - rhs),
            (Value::Constant(constant), Value::Assigned(assigned)) => self
                .scalar_chip()
                .sum_with_coeff_and_const(
                    &mut self.ctx_mut(),
                    &[(-C::Scalar::ONE, assigned)],
                    *constant,
                )
                .map(Value::Assigned)
                .unwrap(),
            (Value::Assigned(assigned), Value::Constant(constant)) => self
                .scalar_chip()
                .sum_with_coeff_and_const(
                    &mut self.ctx_mut(),
                    &[(C::Scalar::ONE, assigned)],
                    -*constant,
                )
                .map(Value::Assigned)
                .unwrap(),
            (Value::Assigned(lhs), Value::Assigned(rhs)) => {
                IntegerInstructions::sub(self.scalar_chip().deref(), &mut self.ctx_mut(), lhs, rhs)
                    .map(Value::Assigned)
                    .unwrap()
            },
        };
        self.scalar(output)
    }

    fn mul(
        self: &Rc<Self>,
        lhs: &Scalar<'a, C, EccChip, ScalarChip>,
        rhs: &Scalar<'a, C, EccChip, ScalarChip>,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        let output = match (lhs.value().deref(), rhs.value().deref()) {
            (Value::Constant(lhs), Value::Constant(rhs)) => Value::Constant(*lhs * rhs),
            (Value::Assigned(assigned), Value::Constant(constant))
            | (Value::Constant(constant), Value::Assigned(assigned)) => self
                .scalar_chip()
                .sum_with_coeff_and_const(
                    &mut self.ctx_mut(),
                    &[(*constant, assigned)],
                    C::Scalar::ZERO,
                )
                .map(Value::Assigned)
                .unwrap(),
            (Value::Assigned(lhs), Value::Assigned(rhs)) => self
                .scalar_chip()
                .sum_products_with_coeff_and_const(
                    &mut self.ctx_mut(),
                    &[(C::Scalar::ONE, lhs, rhs)],
                    C::Scalar::ZERO,
                )
                .map(Value::Assigned)
                .unwrap(),
        };
        self.scalar(output)
    }

    fn neg(
        self: &Rc<Self>,
        scalar: &Scalar<'a, C, EccChip, ScalarChip>,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        let output = match scalar.value().deref() {
            Value::Constant(constant) => Value::Constant(constant.neg()),
            Value::Assigned(assigned) => {
                IntegerInstructions::neg(self.scalar_chip().deref(), &mut self.ctx_mut(), assigned)
                    .map(Value::Assigned)
                    .unwrap()
            },
        };
        self.scalar(output)
    }

    fn invert(
        self: &Rc<Self>,
        scalar: &Scalar<'a, C, EccChip, ScalarChip>,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        let output = match scalar.value().deref() {
            Value::Constant(constant) => Value::Constant(Field::invert(constant).unwrap()),
            Value::Assigned(assigned) => Value::Assigned(
                IntegerInstructions::invert(
                    self.scalar_chip().deref(),
                    &mut self.ctx_mut(),
                    assigned,
                )
                .unwrap(),
            ),
        };
        self.scalar(output)
    }
}

#[cfg(test)]
impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Halo2Loader<'a, C, EccChip, ScalarChip>
{
    fn start_row_metering(self: &Rc<Self>, identifier: &str) {
        use crate::halo2_verifier::loader::halo2::shim::Context;

        self.row_meterings
            .borrow_mut()
            .push((identifier.to_string(), self.ctx().offset()))
    }

    fn end_row_metering(self: &Rc<Self>) {
        use crate::halo2_verifier::loader::halo2::shim::Context;

        let mut row_meterings = self.row_meterings.borrow_mut();
        let (_, row) = row_meterings.last_mut().unwrap();
        *row = self.ctx().offset() - *row;
    }

    pub fn print_row_metering(self: &Rc<Self>) {
        for (identifier, cost) in self.row_meterings.borrow().iter() {
            println!("{identifier}: {cost}");
        }
    }
}

#[derive(Clone, Debug)]
pub enum Value<T, L> {
    Constant(T),
    Assigned(L),
}

impl<T, L> Value<T, L> {
    fn maybe_const(&self) -> Option<T>
    where
        T: Copy,
    {
        match self {
            Value::Constant(constant) => Some(*constant),
            _ => None,
        }
    }

    fn assigned(&self) -> &L {
        match self {
            Value::Assigned(assigned) => assigned,
            _ => unreachable!(),
        }
    }
}

#[derive(Clone)]
pub struct Scalar<
    'a,
    C: CurveAffine,
    EccChip: EccInstructions<'a, C>,
    ScalarChip: IntegerInstructions<'a, C::Base>,
> {
    loader: Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>,
    index: usize,
    value: RefCell<Value<C::Scalar, EccChip::AssignedScalar>>,
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Scalar<'a, C, EccChip, ScalarChip>
{
    pub fn loader(&self) -> &Rc<Halo2Loader<'a, C, EccChip, ScalarChip>> {
        &self.loader
    }

    pub fn into_assigned(self) -> EccChip::AssignedScalar {
        match self.value.into_inner() {
            Value::Constant(constant) => self.loader.assign_const_scalar(constant),
            Value::Assigned(assigned) => assigned,
        }
    }

    pub fn assigned(&self) -> Ref<EccChip::AssignedScalar> {
        if let Some(constant) = self.maybe_const() {
            *self.value.borrow_mut() = Value::Assigned(self.loader.assign_const_scalar(constant))
        }
        Ref::map(self.value.borrow(), Value::assigned)
    }

    fn value(&self) -> Ref<Value<C::Scalar, EccChip::AssignedScalar>> {
        self.value.borrow()
    }

    fn maybe_const(&self) -> Option<C::Scalar> {
        self.value().deref().maybe_const()
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > PartialEq for Scalar<'a, C, EccChip, ScalarChip>
{
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > LoadedScalar<C::Scalar> for Scalar<'a, C, EccChip, ScalarChip>
{
    type Loader = Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>;

    fn loader(&self) -> &Self::Loader {
        &self.loader
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Debug for Scalar<'a, C, EccChip, ScalarChip>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Scalar")
            .field("value", &self.value)
            .finish()
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > FieldOps for Scalar<'a, C, EccChip, ScalarChip>
{
    fn invert(&self) -> Option<Self> {
        Some(self.loader.invert(self))
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Add for Scalar<'a, C, EccChip, ScalarChip>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Halo2Loader::add(&self.loader, &self, &rhs)
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Sub for Scalar<'a, C, EccChip, ScalarChip>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Halo2Loader::sub(&self.loader, &self, &rhs)
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Mul for Scalar<'a, C, EccChip, ScalarChip>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Halo2Loader::mul(&self.loader, &self, &rhs)
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Neg for Scalar<'a, C, EccChip, ScalarChip>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Halo2Loader::neg(&self.loader, &self)
    }
}

impl<
        'a,
        'b,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Add<&'b Self> for Scalar<'a, C, EccChip, ScalarChip>
{
    type Output = Self;

    fn add(self, rhs: &'b Self) -> Self::Output {
        Halo2Loader::add(&self.loader, &self, rhs)
    }
}

impl<
        'a,
        'b,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Sub<&'b Self> for Scalar<'a, C, EccChip, ScalarChip>
{
    type Output = Self;

    fn sub(self, rhs: &'b Self) -> Self::Output {
        Halo2Loader::sub(&self.loader, &self, rhs)
    }
}

impl<
        'a,
        'b,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Mul<&'b Self> for Scalar<'a, C, EccChip, ScalarChip>
{
    type Output = Self;

    fn mul(self, rhs: &'b Self) -> Self::Output {
        Halo2Loader::mul(&self.loader, &self, rhs)
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > AddAssign for Scalar<'a, C, EccChip, ScalarChip>
{
    fn add_assign(&mut self, rhs: Self) {
        *self = Halo2Loader::add(&self.loader, self, &rhs)
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > SubAssign for Scalar<'a, C, EccChip, ScalarChip>
{
    fn sub_assign(&mut self, rhs: Self) {
        *self = Halo2Loader::sub(&self.loader, self, &rhs)
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > MulAssign for Scalar<'a, C, EccChip, ScalarChip>
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = Halo2Loader::mul(&self.loader, self, &rhs)
    }
}

impl<
        'a,
        'b,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > AddAssign<&'b Self> for Scalar<'a, C, EccChip, ScalarChip>
{
    fn add_assign(&mut self, rhs: &'b Self) {
        *self = Halo2Loader::add(&self.loader, self, rhs)
    }
}

impl<
        'a,
        'b,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > SubAssign<&'b Self> for Scalar<'a, C, EccChip, ScalarChip>
{
    fn sub_assign(&mut self, rhs: &'b Self) {
        *self = Halo2Loader::sub(&self.loader, self, rhs)
    }
}

impl<
        'a,
        'b,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > MulAssign<&'b Self> for Scalar<'a, C, EccChip, ScalarChip>
{
    fn mul_assign(&mut self, rhs: &'b Self) {
        *self = Halo2Loader::mul(&self.loader, self, rhs)
    }
}

#[derive(Clone)]
pub struct EcPoint<
    'a,
    C: CurveAffine,
    EccChip: EccInstructions<'a, C>,
    ScalarChip: IntegerInstructions<'a, C::Base>,
> {
    loader: Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>,
    index: usize,
    value: RefCell<Value<C, EccChip::AssignedEcPoint>>,
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > EcPoint<'a, C, EccChip, ScalarChip>
{
    pub fn into_assigned(self) -> EccChip::AssignedEcPoint {
        match self.value.into_inner() {
            Value::Constant(constant) => self.loader.assign_const_ec_point(constant),
            Value::Assigned(assigned) => assigned,
        }
    }

    pub fn assigned(&self) -> Ref<EccChip::AssignedEcPoint> {
        if let Some(constant) = self.maybe_const() {
            *self.value.borrow_mut() = Value::Assigned(self.loader.assign_const_ec_point(constant))
        }
        Ref::map(self.value.borrow(), Value::assigned)
    }

    fn value(&self) -> Ref<Value<C, EccChip::AssignedEcPoint>> {
        self.value.borrow()
    }

    fn maybe_const(&self) -> Option<C> {
        self.value().deref().maybe_const()
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > PartialEq for EcPoint<'a, C, EccChip, ScalarChip>
{
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > LoadedEcPoint<C> for EcPoint<'a, C, EccChip, ScalarChip>
{
    type Loader = Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>;

    fn loader(&self) -> &Self::Loader {
        &self.loader
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Debug for EcPoint<'a, C, EccChip, ScalarChip>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EcPoint")
            .field("index", &self.index)
            .field("value", &self.value)
            .finish()
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > ScalarLoader<C::Scalar> for Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>
{
    type LoadedScalar = Scalar<'a, C, EccChip, ScalarChip>;

    fn load_const(&self, value: &C::Scalar) -> Scalar<'a, C, EccChip, ScalarChip> {
        self.scalar(Value::Constant(*value))
    }

    fn assert_eq(
        &self,
        annotation: &str,
        lhs: &Scalar<'a, C, EccChip, ScalarChip>,
        rhs: &Scalar<'a, C, EccChip, ScalarChip>,
    ) -> Result<(), crate::halo2_verifier::Error> {
        self.scalar_chip()
            .assert_equal(&mut self.ctx_mut(), &lhs.assigned(), &rhs.assigned())
            .map_err(|_| crate::halo2_verifier::Error::AssertionFailure(annotation.to_string()))
    }

    fn sum_with_coeff_and_const(
        &self,
        values: &[(C::Scalar, &Scalar<'a, C, EccChip, ScalarChip>)],
        constant: C::Scalar,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        let values = values
            .iter()
            .map(|(coeff, value)| (*coeff, value.assigned()))
            .collect_vec();
        self.scalar(Value::Assigned(
            self.scalar_chip()
                .sum_with_coeff_and_const(&mut self.ctx_mut(), &values, constant)
                .unwrap(),
        ))
    }

    fn sum_products_with_coeff_and_const(
        &self,
        values: &[(
            C::Scalar,
            &Scalar<'a, C, EccChip, ScalarChip>,
            &Scalar<'a, C, EccChip, ScalarChip>,
        )],
        constant: C::Scalar,
    ) -> Scalar<'a, C, EccChip, ScalarChip> {
        let values = values
            .iter()
            .map(|(coeff, lhs, rhs)| (*coeff, lhs.assigned(), rhs.assigned()))
            .collect_vec();
        self.scalar(Value::Assigned(
            self.scalar_chip()
                .sum_products_with_coeff_and_const(&mut self.ctx_mut(), &values, constant)
                .unwrap(),
        ))
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > EcPointLoader<C> for Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>
{
    type LoadedEcPoint = EcPoint<'a, C, EccChip, ScalarChip>;

    fn ec_point_load_const(&self, ec_point: &C) -> EcPoint<'a, C, EccChip, ScalarChip> {
        self.ec_point(Value::Constant(*ec_point))
    }

    fn ec_point_assert_eq(
        &self,
        annotation: &str,
        lhs: &EcPoint<'a, C, EccChip, ScalarChip>,
        rhs: &EcPoint<'a, C, EccChip, ScalarChip>,
    ) -> Result<(), crate::halo2_verifier::Error> {
        if let (Value::Constant(lhs), Value::Constant(rhs)) =
            (lhs.value().deref(), rhs.value().deref())
        {
            assert_eq!(lhs, rhs);
            Ok(())
        } else {
            let lhs = lhs.assigned();
            let rhs = rhs.assigned();
            self.ecc_chip()
                .assert_equal(&mut self.ctx_mut(), lhs.deref(), rhs.deref())
                .map_err(|_| crate::halo2_verifier::Error::AssertionFailure(annotation.to_string()))
        }
    }

    fn multi_scalar_multiplication(
        pairs: &[(
            &<Self as ScalarLoader<C::Scalar>>::LoadedScalar,
            &EcPoint<'a, C, EccChip, ScalarChip>,
        )],
    ) -> EcPoint<'a, C, EccChip, ScalarChip> {
        let loader = &pairs[0].0.loader;

        let (constant, fixed_base, variable_base_non_scaled, variable_base_scaled) =
            pairs.iter().cloned().fold(
                (C::identity(), Vec::new(), Vec::new(), Vec::new()),
                |(
                    mut constant,
                    mut fixed_base,
                    mut variable_base_non_scaled,
                    mut variable_base_scaled,
                ),
                 (scalar, base)| {
                    match (scalar.value().deref(), base.value().deref()) {
                        (Value::Constant(scalar), Value::Constant(base)) => {
                            constant = (*base * scalar + constant).into()
                        },
                        (Value::Assigned(_), Value::Constant(base)) => {
                            fixed_base.push((scalar, *base))
                        },
                        (Value::Constant(scalar), Value::Assigned(_))
                            if scalar.eq(&C::Scalar::ONE) =>
                        {
                            variable_base_non_scaled.push(base);
                        },
                        _ => variable_base_scaled.push((scalar, base)),
                    };
                    (
                        constant,
                        fixed_base,
                        variable_base_non_scaled,
                        variable_base_scaled,
                    )
                },
            );

        let fixed_base_msm = (!fixed_base.is_empty())
            .then(|| {
                let fixed_base = fixed_base
                    .into_iter()
                    .map(|(scalar, base)| (scalar.assigned(), base))
                    .collect_vec();
                loader
                    .ecc_chip
                    .borrow_mut()
                    .fixed_base_msm(&mut loader.ctx_mut(), &fixed_base)
                    .unwrap()
            })
            .map(RefCell::new);
        let variable_base_msm = (!variable_base_scaled.is_empty())
            .then(|| {
                let variable_base_scaled = variable_base_scaled
                    .into_iter()
                    .map(|(scalar, base)| (scalar.assigned(), base.assigned()))
                    .collect_vec();
                loader
                    .ecc_chip
                    .borrow_mut()
                    .variable_base_msm(&mut loader.ctx_mut(), &variable_base_scaled)
                    .unwrap()
            })
            .map(RefCell::new);
        let output = loader
            .ecc_chip()
            .sum_with_const(
                &mut loader.ctx_mut(),
                &variable_base_non_scaled
                    .into_iter()
                    .map(EcPoint::assigned)
                    .chain(fixed_base_msm.as_ref().map(RefCell::borrow))
                    .chain(variable_base_msm.as_ref().map(RefCell::borrow))
                    .collect_vec(),
                constant,
            )
            .unwrap();

        loader.ec_point_from_assigned(output)
    }
}

impl<
        'a,
        C: CurveAffine,
        EccChip: EccInstructions<'a, C>,
        ScalarChip: IntegerInstructions<'a, C::Base>,
    > Loader<C> for Rc<Halo2Loader<'a, C, EccChip, ScalarChip>>
{
    #[cfg(test)]
    fn start_cost_metering(&self, identifier: &str) {
        self.start_row_metering(identifier)
    }

    #[cfg(test)]
    fn end_cost_metering(&self) {
        self.end_row_metering()
    }
}
