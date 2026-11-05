use core::ops::{Add, Mul, Neg, Sub};
use rand::Rng;
use std::{
    fmt::Debug,
    ops::{AddAssign, MulAssign, SubAssign},
};

pub trait FiniteField:
    Sized
    + Default
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + Copy
    + Clone
    + PartialEq
{
    fn zero() -> Self;

    fn one() -> Self;

    fn random<R: Rng + ?Sized>(rng: &mut R) -> Self;

    fn expand_from_seed(seed: &[u8], count: usize) -> Vec<Self>;

    fn pow(self, mut exp: u64) -> Self {
        let mut base = self;
        let mut acc = Self::one();

        while exp > 0 {
            if (exp & 1) == 1 {
                acc *= base;
            }
            base *= base;
            exp >>= 1;
        }
        acc
    }

    fn inv(&self) -> Self {
        /* To implement...!!! */
        Self::one()
    }
}

pub trait ExtensionConfig<const DEGREE: usize>: Copy + Default + PartialEq + Debug {
    type Field: FiniteField;

    const REDUCE_COEFFS: [Self::Field; DEGREE];
}
