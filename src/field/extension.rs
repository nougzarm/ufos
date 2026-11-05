//! # General Finite Fields Operations
//!
//! This module provides the [`Fq`] structure, which allows for operations on finite fields
//! by implementing the [`FiniteField`] trait.

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

/// Finite field (As extension of prime field)
///
/// Implements the finite field $\mathbb{F}_q$ containing $q = p^e$ elements, where $p$ is a prime number.
/// Such a finite field can be described as an extension of prime finite field.
///
/// # Type Parameters
///
/// * `P` - Prime field characteristic ($\mathbb{F}_p$).
/// * `E` - Extension degree ($\mathbb{F}_{p^e} = \mathbb{F}_q$).
/// * `C` - Field extension parameters.
///
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
    /// Creates a new extension field element from its canonical basis coefficients.
    ///
    /// # Arguments
    ///
    /// * `coeffs` - Array of $e$ coefficients in $\mathbb{F}_p$ in ascending degree order $(c_0, \dots, c_{e-1})$.
    pub const fn new(coeffs: [Fp<P>; E]) -> Self {
        Self {
            coeffs,
            _config: PhantomData,
        }
    }

    /// Returns the additive identity $0 \in \mathbb{F}_q$.
    pub fn zero() -> Self {
        Self::new([<Fp<P> as FiniteField>::zero(); E])
    }

    /// Returns the canonical generator $\xi$ of the extension.
    ///
    /// Represented as the polynomial $0 + 1 \cdot \xi^1 + 0 \cdot \xi^2 + \dots + 0 \cdot \xi^{e-1}$.
    pub fn generator() -> Self {
        let mut xi = Self::zero();
        xi.coeffs[1] = Fp::<P>::one();
        xi
    }

    /// Returns the canonical basis $(1, \xi, \xi^2, \dots, \xi^{e-1})$ of $\mathbb{F}_q$ over $\mathbb{F}_p$.
    ///
    /// The $i$-th element of the returned array has coefficient $1$ at index $i$ and $0$ elsewhere.
    pub fn base() -> [Self; E] {
        let res: [Self; E] = core::array::from_fn(|i| {
            let mut coeffs = [Fp::<P>::zero(); E];
            coeffs[i] = Fp::<P>::one();
            Self::new(coeffs)
        });
        res
    }

    /// Computes the $p$-th power of the canonical generator $\xi$.
    ///
    /// Used as the primary building block for Horner-based evaluation of the Frobenius endomorphism.
    pub fn xi_pow_p() -> Self {
        let xi = Self::generator();
        xi.pow(P as u64)
    }

    /// Evaluates the Frobenius automorphism $\phi(\alpha) = \alpha^p$.
    ///
    /// Exploits the $\mathbb{F}_p$-linearity of the Frobenius map (since $c_j^p = c_j$ for all $c_j \in \mathbb{F}_p$):
    ///
    /// $$\phi\left(\sum_{j=0}^{e-1} c_j \xi^j\right) = \sum_{j=0}^{e-1} c_j (\xi^p)^j$$
    ///
    /// Evaluated efficiently using Horner's rule with $\xi^p$.
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

    /// Computes the iterated Frobenius power $\phi^a(\alpha) = \alpha^{p^a}$.
    ///
    /// Since the Galois group $\mathrm{Gal}(\mathbb{F}_{p^e}/\mathbb{F}_p)$ is cyclic of order $e$,
    /// the exponent $a$ is automatically reduced modulo $e$:
    ///
    /// $$\alpha^{[a]} = \alpha^{p^{a \pmod e}}$$
    ///
    /// # Arguments
    ///
    /// * `a` - Frobenius exponent.
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
    /// Embeds a base field element $c \in \mathbb{F}_p$ into the extension field $\mathbb{F}_q$ as $c \cdot \xi^0$.
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

    /// Computes component-wise addition in $\mathbb{F}_{p^e}$.
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

    /// Computes the additive inverse by negating each coefficient in $\mathbb{F}_p$.
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

    /// Computes component-wise subtraction in $\mathbb{F}_{p^e}$.
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

    /// Multiplies two extension field elements modulo the irreducible polynomial $i(\xi)$.
    ///
    /// Computes the unreduced polynomial product of degree at most $2e - 2$, split into
    /// low-degree (`low`) and high-degree (`high`) parts, followed by iterative reduction
    /// using `C::REDUCE_COEFFS`.
    fn mul(self, rhs: Self) -> Self::Output {
        let mut low = [<Fp<P> as FiniteField>::zero(); E];
        let mut high = [<Fp<P> as FiniteField>::zero(); E];

        for i in 0..E {
            let split = E - i;

            for j in 0..split {
                low[i + j] += self.coeffs[i] * rhs.coeffs[j];
            }

            for j in split..E {
                high[i + j - E] += self.coeffs[i] * rhs.coeffs[j];
            }
        }

        for d in (0..E - 1).rev() {
            let c = high[d];
            let split = E - d;

            for j in 0..split {
                low[d + j] += c * C::REDUCE_COEFFS[j];
            }

            for j in split..E {
                high[d + j - E] += c * C::REDUCE_COEFFS[j];
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
    /// Returns the additive identity $0 \in \mathbb{F}_q$.
    fn zero() -> Self {
        Self::zero()
    }

    /// Returns the multiplicative identity $1 \in \mathbb{F}_q$.
    fn one() -> Self {
        Self::from(<Fp<P> as FiniteField>::one())
    }

    /// Generates a uniformly random field element by sampling $e$ random coordinates in $\mathbb{F}_p$.
    ///
    /// # Arguments
    ///
    /// * `rng` - Pseudo-random generator reference.
    fn random<R: Rng + ?Sized>(rng: &mut R) -> Self {
        let coeffs: [Fp<P>; E] = core::array::from_fn(|_| Fp::random(rng));
        Self::new(coeffs)
    }

    /// Deterministically expands a byte seed into a sequence of extension field elements using `Shake256`.
    ///
    /// Expands `count * E` base field elements from `seed` and groups them into $e$-tuples.
    ///
    /// # Arguments
    ///
    /// * `seed` - The input entropy byte slice.
    /// * `count` - The number of [`Fq`] elements to generate.
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
    use super::*;

    const P: u8 = 61;
    type F = crate::params::F31_61;

    #[test]
    fn fq_field_laws() {
        for _ in 0..P {
            let a = F::random(&mut rand_core::OsRng);
            assert_eq!(a + F::zero(), a);
            assert_eq!(a * F::one(), a);
            assert_eq!(a + (-a), F::zero());

            for _ in 0..P {
                let b = F::random(&mut rand_core::OsRng);
                assert_eq!(a + b, b + a);
                assert_eq!(a * b, b * a);

                for _ in 0..P {
                    let c = F::random(&mut rand_core::OsRng);
                    assert_eq!(a * (b + c), a * b + a * c);
                }
            }
        }
    }
}
