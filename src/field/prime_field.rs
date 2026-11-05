//! # Prime Fields Operations
//!
//! This module provides the [`Fp`] structure, which allows for operations on basic prime finite fields
//! that implement the [`FiniteField`] trait.

use core::ops::{Add, Mul, MulAssign, Sub};
use rand::Rng;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use std::ops::{AddAssign, Neg, SubAssign};

use crate::field::traits::FiniteField;

/// Prime finite field
///
/// Implements the finite field $\mathbb{F}_p$ containing $p$ elements, where $p$ is a prime number.
/// Such a field can be described as the ring $\mathbb{Z}/p\mathbb{Z}$.
///
/// # Type Parameters
///
/// * `P` - Field prime cardinal $p$.
///
/// # Examples
///
/// ```rust
/// use ufos::field::{
///     traits::FiniteField,
///     prime_field::Fp
/// };
///
/// type F = Fp<31>;
/// let one = F::one();
/// let x = F::from(3) * one;
/// ```
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub struct Fp<const P: u8>(u8);

impl<const P: u8> Fp<P> {
    const CHECK_P: () = assert!(P > 1 && P < (1 << 7), "P must be in 2..2^7");

    /// Creates a field element from a raw `u8` integer, reducing it modulo $p$.
    ///
    /// # Arguments
    ///
    /// * `val` - An unsigned 8-bit integer.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use ufos::field::prime_field::Fp;
    ///
    /// type F = Fp<31>;
    /// let a = F::new(35);
    /// assert_eq!(a.value(), 4);
    /// ```
    pub const fn new(val: u8) -> Self {
        let _ = Self::CHECK_P;
        Self(val % P)
    }

    /// Instantiates a field element assuming `val` is already in canonical form ($0 \le \mathrm{val} < p$).
    ///
    /// This bypasses modular division and is intended for performance-critical inner loops
    /// when the input is guaranteed to be reduced.
    ///
    /// # Arguments
    ///
    /// * `val` - A pre-reduced integer strictly less than $p$.
    pub const fn from_canonical(val: u8) -> Self {
        Self(val)
    }

    /// Returns the canonical integer representation of the field element in $[0, p - 1]$.
    pub const fn value(self) -> u8 {
        self.0
    }
}

impl<const P: u8> From<u8> for Fp<P> {
    fn from(value: u8) -> Self {
        Self(value % P)
    }
}

impl<const P: u8> Add for &Fp<P> {
    type Output = Fp<P>;

    /// Computes $(a + b) \pmod p$ without integer division.
    fn add(self, rhs: Self) -> Self::Output {
        let sum = self.0 + rhs.0;
        let reduced = if sum >= P { sum - P } else { sum };
        Fp::<P>::from_canonical(reduced)
    }
}

impl<const P: u8> Add for Fp<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl<const P: u8> AddAssign for Fp<P> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const P: u8> Neg for &Fp<P> {
    type Output = Fp<P>;

    /// Computes the additive inverse $-a \pmod p$.
    ///
    /// Maps $0 \mapsto 0$ and $a \mapsto p - a$ for $a > 0$.
    fn neg(self) -> Self::Output {
        Fp::<P>::from(P - self.value())
    }
}

impl<const P: u8> Neg for Fp<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<const P: u8> Sub for &Fp<P> {
    type Output = Fp<P>;

    /// Computes $(a - b) \pmod p$ with conditional borrow correction.
    fn sub(self, rhs: Self) -> Self::Output {
        let (diff, borrow) = self.0.overflowing_sub(rhs.0);
        let corr = if borrow { P } else { 0 };
        Fp::<P>(diff.wrapping_add(corr))
    }
}

impl<const P: u8> Sub for Fp<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl<const P: u8> SubAssign for Fp<P> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const P: u8> Mul for &Fp<P> {
    type Output = Fp<P>;

    /// Computes $(a * b) \pmod p$ using a widening 16-bit intermediate product.
    fn mul(self, rhs: Self) -> Self::Output {
        let prod = (self.0 as u16) * (rhs.0 as u16);

        Fp::<P>((prod % (P as u16)) as u8)
    }
}

impl<const P: u8> Mul for Fp<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl<const P: u8> MulAssign for Fp<P> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const P: u8> FiniteField for Fp<P> {
    /// Returns the additive identity $0 \in \mathbb{F}_p$.
    #[inline]
    fn zero() -> Self {
        Self::from_canonical(0u8)
    }

    /// Returns the multiplicative identity $1 \in \mathbb{F}_p$.
    #[inline]
    fn one() -> Self {
        Self::from_canonical(1u8)
    }

    /// Samples a uniformly distributed random element from $\mathbb{F}_p$.
    ///
    /// # Arguments
    ///
    /// * `rng` - Pseudo-random generator reference.
    fn random<R: Rng + ?Sized>(rng: &mut R) -> Self {
        let val = rng.gen_range(0..P);
        Self::from_canonical(val)
    }

    /// Deterministically expands a byte seed into a sequence of field elements using `Shake256`.
    ///
    /// # Arguments
    ///
    /// * `seed` - The input entropy byte slice.
    /// * `count` - The number of field elements to generate.
    fn expand_from_seed(seed: &[u8], count: usize) -> Vec<Self> {
        let mut hasher = Shake256::default();
        hasher.update(seed);
        let mut reader = hasher.finalize_xof();

        let mut list = Vec::with_capacity(count);
        let mut buffer = [0u8; 1];

        while list.len() < count {
            reader.read(&mut buffer);
            let val = buffer[0];
            list.push(Self::from(val));
        }

        list
    }

    /// Computes the multiplicative inverse $a^{-1} \pmod p$ using Fermat's Little Theorem:
    ///
    /// $$a^{-1} \equiv a^{p - 2} \pmod p$$
    ///
    /// # Panics
    ///
    /// Panics if attempting to invert zero ($a = 0$).
    fn inv(&self) -> Self {
        assert!(self.0 != 0);

        let p_val = P as u64;
        self.pow(p_val - 2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: u8 = 61;
    type F = Fp<P>;
    type F2 = Fp<2>;
    const MERSENNE_7: u8 = (1 << 7) - 1;
    type FMersenne = Fp<MERSENNE_7>;

    #[test]
    fn fp_field_laws() {
        for x in 0..P {
            let a = F::new(x);
            assert_eq!(a + F::zero(), a);
            assert_eq!(a * F::one(), a);
            assert_eq!(a + (-a), F::zero());

            if a != F::zero() {
                assert_eq!(a * a.inv(), F::one());
            }

            for y in 0..P {
                let b = F::new(y);
                assert_eq!(a + b, b + a);
                assert_eq!(a * b, b * a);

                for z in 0..P {
                    let c = F::new(z);
                    assert_eq!(a * (b + c), a * b + a * c);
                }
            }
        }
    }

    #[test]
    fn test_fermat_little_theorem_f31() {
        for a_val in 1..P {
            let a = F::new(a_val);
            let mut pow = F::one();
            for _ in 0..(P - 1) {
                pow *= a;
            }
            assert_eq!(pow, F::one());
        }
    }

    #[test]
    fn test_boundary_max_prime_near_2_pow_31() {
        let max_elem = FMersenne::new(MERSENNE_7 - 1);
        let one = FMersenne::one();
        let zero = FMersenne::zero();

        // (P - 1) + 1 == 0 mod P
        assert_eq!(max_elem + one, zero);

        // (P - 1) + (P - 1) == P - 2 mod P
        assert_eq!((max_elem + max_elem).value(), MERSENNE_7 - 2);

        // 0 - 1 == P - 1 mod P
        assert_eq!(zero - one, max_elem);

        // 0 - (P - 1) == 1 mod P
        assert_eq!(zero - max_elem, one);

        // (-1) * (-1) == 1 mod P
        assert_eq!(max_elem * max_elem, one);
    }

    #[test]
    fn test_f2_edge_case() {
        let z = F2::zero();
        let o = F2::one();

        assert_eq!(o + o, z);
        assert_eq!(z - o, o);
        assert_eq!(o * o, o);
        assert_eq!(o * z, z);
    }
}
