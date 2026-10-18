use core::ops::{Add, Mul, MulAssign, Sub};

use crate::field::traits::FiniteField;

/* P should be < 2^7 for the Add trait implementation */
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub struct Fp<const P: u8>(u8);

impl<const P: u8> Fp<P> {
    const CHECK_P: () = assert!(P > 1 && P < (1 << 7), "P must be in 2..2^7");

    pub const fn new(val: u8) -> Self {
        let _ = Self::CHECK_P;
        Self(val % P)
    }

    pub const fn from_canonical(val: u8) -> Self {
        Self(val)
    }

    pub const fn value(self) -> u8 {
        self.0
    }
}

impl<const P: u8> From<u8> for Fp<P> {
    fn from(value: u8) -> Self {
        Self(value % P)
    }
}

impl<const P: u8> Add for Fp<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let sum = self.0 + rhs.0;
        let reduced = if sum >= P { sum - P } else { sum };
        Self(reduced)
    }
}

impl<const P: u8> Sub for Fp<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let (diff, borrow) = self.0.overflowing_sub(rhs.0);
        let corr = if borrow { P } else { 0 };
        Self(diff.wrapping_add(corr))
    }
}

impl<const P: u8> Mul for Fp<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let prod = (self.0 as u16) * (rhs.0 as u16);

        Self((prod % (P as u16)) as u8)
    }
}

impl<const P: u8> MulAssign for Fp<P> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const P: u8> FiniteField for Fp<P> {
    fn zero() -> Self {
        Self::new(0u8)
    }

    fn one() -> Self {
        Self::new(1u8)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type F31 = Fp<31>;
    type F2 = Fp<2>;
    const MERSENNE_7: u8 = (1 << 7) - 1;
    type FMersenne = Fp<MERSENNE_7>;

    #[test]
    fn test_canonical_reduction_on_new() {
        assert_eq!(F31::new(0).value(), 0);
        assert_eq!(F31::new(1).value(), 1);
        assert_eq!(F31::new(32).value(), 1);
        assert_eq!(F31::new(u8::MAX).value(), u8::MAX % 31);
    }

    #[test]
    fn test_identities_and_default() {
        let zero = F31::zero();
        let one = F31::one();
        assert_eq!(zero, F31::default());
        assert_eq!(one.value(), 1);

        for a_val in 0..31 {
            let a = F31::new(a_val);
            assert_eq!(a + zero, a);
            assert_eq!(zero + a, a);
            assert_eq!(a - zero, a);
            assert_eq!(a - a, zero);
            assert_eq!(a * one, a);
            assert_eq!(one * a, a);
            assert_eq!(a * zero, zero);
        }
    }

    #[test]
    fn test_exhaustive_ring_axioms_f31() {
        for a_val in 0..31 {
            let a = F31::from_canonical(a_val);
            for b_val in 0..31 {
                let b = F31::from_canonical(b_val);

                // Vérification contre l'arithmétique u16 de référence
                assert_eq!((a + b).value(), (a_val + b_val) % 31);
                assert_eq!((a - b).value(), (a_val + 31 - b_val) % 31);
                assert_eq!((a * b).value(), (a_val * b_val) % 31);

                // Commutativity
                assert_eq!(a + b, b + a);
                assert_eq!(a * b, b * a);

                // Cohérence Add / Sub
                assert_eq!((a + b) - b, a);
                assert_eq!((a - b) + b, a);

                // MulAssign
                let mut acc = a;
                acc *= b;
                assert_eq!(acc, a * b);

                for c_val in 0..31 {
                    let c = F31::from_canonical(c_val);
                    // Associativité
                    assert_eq!((a + b) + c, a + (b + c));
                    assert_eq!((a * b) * c, a * (b * c));
                    // Distributivité
                    assert_eq!(a * (b + c), (a * b) + (a * c));
                    assert_eq!(a * (b - c), (a * b) - (a * c));
                }
            }
        }
    }

    #[test]
    fn test_fermat_little_theorem_f31() {
        // Pour tout a != 0 dans F_p, a^(p-1) == 1
        for a_val in 1..31 {
            let a = F31::new(a_val);
            let mut pow = F31::one();
            for _ in 0..30 {
                pow *= a;
            }
            assert_eq!(
                pow,
                F31::one(),
                "Échec du petit théorème de Fermat pour a = {}",
                a_val
            );
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
