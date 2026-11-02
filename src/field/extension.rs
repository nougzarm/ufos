use core::ops::{Add, Mul, MulAssign, Sub};
use rand::Rng;
use std::{
    marker::PhantomData,
    ops::{AddAssign, Neg, SubAssign},
};

use crate::field::{
    prime_field::Fp,
    traits::{ExtensionConfig, FiniteField},
};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Fq<const P: u8, const E: usize, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    pub coeffs: [C::Field; E],
    _config: PhantomData<C>,
}

impl<const P: u8, const E: usize, C> Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
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

    pub fn generator() -> Self {
        let mut xi = Self::zero();
        xi.coeffs[1] = Fp::<P>::one();
        xi
    }

    pub fn base() -> [Self; E] {
        let res: [Self; E] = core::array::from_fn(|i| {
            let mut coeffs = [Fp::<P>::zero(); E];
            coeffs[i] = Fp::<P>::one();
            Self::new(coeffs)
        });

        res
    }

    /* Returns xi^p where xi is the canonical generator */
    pub fn xi_pow_p() -> Self {
        let xi = Self::generator();
        xi.pow(P as u64)
    }

    pub fn frobenius(&self) -> Self {
        let xi_p = Self::xi_pow_p();
        let mut res = Self::zero();

        for i in (0..E).rev() {
            if i == E - 1 {
                res = Self::from(self.coeffs[i]);
            } else {
                res = res * xi_p + Self::from(self.coeffs[i]);
            }
        }
        res
    }

    pub fn frobenius_pow(&self, a: usize) -> Self {
        let mut res = *self;
        let shifts = a % E;
        for _ in 0..shifts {
            res = res.frobenius();
        }
        res
    }
}

impl<const P: u8, const E: usize, C> Default for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    fn default() -> Self {
        Self::zero()
    }
}

impl<const P: u8, const E: usize, C> From<Fp<P>> for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    fn from(value: Fp<P>) -> Self {
        let mut result = Self::zero();
        result.coeffs[0] = value;
        result
    }
}

impl<const P: u8, const E: usize, C> Add for &Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    type Output = Fq<P, E, C>;

    fn add(self, rhs: Self) -> Self::Output {
        Fq::<P, E, C>::new(core::array::from_fn(|i| self.coeffs[i] + rhs.coeffs[i]))
    }
}

impl<const P: u8, const E: usize, C> Add for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl<const P: u8, const E: usize, C> AddAssign for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    fn add_assign(&mut self, rhs: Self) {
        for (i, item) in self.coeffs.iter_mut().enumerate() {
            *item += rhs.coeffs[i];
        }
    }
}

impl<const P: u8, const E: usize, C> Neg for &Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    type Output = Fq<P, E, C>;

    fn neg(self) -> Self::Output {
        Fq::<P, E, C>::new(core::array::from_fn(|i| -&self.coeffs[i]))
    }
}

impl<const P: u8, const E: usize, C> Neg for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<const P: u8, const E: usize, C> Sub for &Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    type Output = Fq<P, E, C>;

    fn sub(self, rhs: Self) -> Self::Output {
        Fq::<P, E, C>::new(core::array::from_fn(|i| self.coeffs[i] - rhs.coeffs[i]))
    }
}

impl<const P: u8, const E: usize, C> Sub for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl<const P: u8, const E: usize, C> SubAssign for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    fn sub_assign(&mut self, rhs: Self) {
        for (i, item) in self.coeffs.iter_mut().enumerate() {
            *item -= rhs.coeffs[i];
        }
    }
}

impl<const P: u8, const E: usize, C> Mul for &Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    type Output = Fq<P, E, C>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut low = [<Fp<P> as FiniteField>::zero(); E];
        let mut high = [<Fp<P> as FiniteField>::zero(); E];

        for i in 0..E {
            for j in 0..E {
                let prod = self.coeffs[i] * rhs.coeffs[j];
                if i + j < E {
                    low[i + j] += prod;
                } else {
                    high[i + j - E] += prod;
                }
            }
        }

        for d in (0..E - 1).rev() {
            let c = high[d];
            for j in 0..E {
                let term = c * C::REDUCE_COEFFS[j];
                if d + j < E {
                    low[d + j] += term;
                } else {
                    high[d + j - E] += term;
                }
            }
        }

        Fq::<P, E, C>::new(low)
    }
}

impl<const P: u8, const E: usize, C> Mul for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl<const P: u8, const E: usize, C> MulAssign for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const P: u8, const E: usize, C> FiniteField for Fq<P, E, C>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    fn zero() -> Self {
        Self::zero()
    }

    fn one() -> Self {
        Self::from(<Fp<P> as FiniteField>::one())
    }

    fn random<R: Rng + ?Sized>(rng: &mut R) -> Self {
        let coeffs: [Fp<P>; E] = core::array::from_fn(|_| Fp::random(rng));
        Self::new(coeffs)
    }

    fn expand_from_seed(seed: &[u8], count: usize) -> Vec<Self> {
        let coeffs_list = <Fp<P> as FiniteField>::expand_from_seed(seed, count * E);
        let mut list = Vec::<Self>::with_capacity(count);

        for i in 0..count {
            let array: [Fp<P>; E] = core::array::from_fn(|j| coeffs_list[i * E + j]);
            list.push(Self::new(array));
        }

        list
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
        let r = F31_61::random(&mut rand::rngs::OsRng);

        assert_eq!(r, r * one);
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
