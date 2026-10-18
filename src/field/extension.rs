use core::ops::{Add, Mul, MulAssign, Sub};
use std::marker::PhantomData;

use crate::field::{
    prime_field::Fp,
    traits::{ExtensionConfig, FiniteField},
};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Fq<const P: u8, const E: usize, C>
where
    C: ExtensionConfig<Fp<P>, E>,
{
    pub coeffs: [Fp<P>; E],
    _config: PhantomData<C>,
}

impl<const P: u8, const E: usize, C> Fq<P, E, C>
where
    C: ExtensionConfig<Fp<P>, E>,
{
    pub const fn new(coeffs: [Fp<P>; E]) -> Self {
        Self {
            coeffs,
            _config: PhantomData,
        }
    }

    pub fn zero() -> Self {
        Self::new([<Fp<P> as FiniteField>::zero(); E])
    }
}

impl<const P: u8, const E: usize, C> From<Fp<P>> for Fq<P, E, C>
where
    C: ExtensionConfig<Fp<P>, E>,
{
    fn from(value: Fp<P>) -> Self {
        let mut result = Self::zero();
        result.coeffs[0] = value;
        result
    }
}

impl<const P: u8, const E: usize, C> Add for Fq<P, E, C>
where
    C: ExtensionConfig<Fp<P>, E>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(core::array::from_fn(|i| self.coeffs[i] + rhs.coeffs[i]))
    }
}

impl<const P: u8, const E: usize, C> Sub for Fq<P, E, C>
where
    C: ExtensionConfig<Fp<P>, E>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(core::array::from_fn(|i| self.coeffs[i] - rhs.coeffs[i]))
    }
}

impl<const P: u8, const E: usize, C> Mul for Fq<P, E, C>
where
    C: ExtensionConfig<Fp<P>, E>,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut low = [<Fp<P> as FiniteField>::zero(); E];
        let mut high = [<Fp<P> as FiniteField>::zero(); E];

        for i in 0..E {
            for j in 0..E {
                let prod = self.coeffs[i] * rhs.coeffs[j];
                if i + j < E {
                    low[i + j] = low[i + j] + prod;
                } else {
                    high[i + j - E] = high[i + j - E] + prod;
                }
            }
        }

        for d in (0..E - 1).rev() {
            let c = high[d];
            for j in 0..E {
                let term = c * C::REDUCE_COEFFS[j];
                if d + j < E {
                    low[d + j] = low[d + j] + term;
                } else {
                    high[d + j - E] = high[d + j - E] + term;
                }
            }
        }

        Self::new(low)
    }
}

impl<const P: u8, const E: usize, C> MulAssign for Fq<P, E, C>
where
    C: ExtensionConfig<Fp<P>, E>,
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const P: u8, const E: usize, C> FiniteField for Fq<P, E, C>
where
    C: ExtensionConfig<Fp<P>, E>,
{
    fn zero() -> Self {
        Self::zero()
    }

    fn one() -> Self {
        Self::from(<Fp<P> as FiniteField>::one())
    }
}

#[cfg(test)]
mod ufs_tests {
    use std::assert_ne;

    use super::*;
    use crate::field::prime_field::Fp;
    use crate::params::F31_61;

    #[test]
    fn test_ufo_level1_from_params() {
        let zero = F31_61::zero();
        let one = F31_61::one();

        assert_eq!(zero * one, zero);
        assert_eq!(one * one, one);
        assert_eq!(zero + one, one);
        assert_ne!(zero + one, zero);
    }

    #[test]
    fn test_ufo_level1_multiplication() {
        let mut a = F31_61::zero();
        a.coeffs[0] = Fp::new(5);
        a.coeffs[2] = Fp::new(12);

        let one = F31_61::one();
        assert_eq!(a * one, a);
    }
}
