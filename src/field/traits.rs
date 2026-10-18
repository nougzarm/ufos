use core::ops::{Add, Mul, Sub};
use std::ops::MulAssign;

pub trait FiniteField: Sized + Add + Sub + Mul + MulAssign + Copy + Clone {
    fn zero() -> Self;

    fn one() -> Self;

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
}

pub trait PrimeField: FiniteField {
    fn as_integer(&self) -> [u8];
}

pub trait ExtensionConfig<F: FiniteField, const E: usize>: Copy {
    const REDUCE_COEFFS: [F; E];
}
