use ufos::{
    field::{
        extension::Fq,
        prime_field::Fp,
        traits::FiniteField,
    },
    params::{F31_61, Level1Config},
};

type F = F31_61;
type Base = Fp<31>;

fn element(coeffs: [u8; 61]) -> F {
    F::new(core::array::from_fn(|i| Base::new(coeffs[i])))
}

#[test]
fn extension_addition_and_negation() {
    let a = element(core::array::from_fn(|i| (i % 31) as u8));
    let b = element(core::array::from_fn(|i| ((3 * i + 7) % 31) as u8));

    assert_eq!(a + F::zero(), a);
    assert_eq!(a + (-a), F::zero());
    assert_eq!(a + b, b + a);
    assert_eq!((a + b) - b, a);
}

#[test]
fn extension_multiplication_identity_and_zero() {
    let a = element(core::array::from_fn(|i| ((5 * i + 2) % 31) as u8));

    assert_eq!(a * F::one(), a);
    assert_eq!(F::one() * a, a);
    assert_eq!(a * F::zero(), F::zero());
    assert_eq!(F::zero() * a, F::zero());
}

#[test]
fn extension_multiplication_is_commutative() {
    let a = element(core::array::from_fn(|i| ((7 * i + 1) % 31) as u8));
    let b = element(core::array::from_fn(|i| ((11 * i + 4) % 31) as u8));

    assert_eq!(a * b, b * a);
}

#[test]
fn extension_multiplication_is_associative() {
    let a = element(core::array::from_fn(|i| ((i * 3 + 1) % 31) as u8));
    let b = element(core::array::from_fn(|i| ((i * 5 + 2) % 31) as u8));
    let c = element(core::array::from_fn(|i| ((i * 9 + 6) % 31) as u8));

    assert_eq!((a * b) * c, a * (b * c));
}

#[test]
fn extension_multiplication_distributes_over_addition() {
    let a = element(core::array::from_fn(|i| ((i * 2 + 3) % 31) as u8));
    let b = element(core::array::from_fn(|i| ((i * 7 + 1) % 31) as u8));
    let c = element(core::array::from_fn(|i| ((i * 13 + 8) % 31) as u8));

    assert_eq!(a * (b + c), a * b + a * c);
    assert_eq!((a + b) * c, a * c + b * c);
}

#[test]
fn extension_base_field_embedding_multiplies_correctly() {
    let x = Base::new(17);
    let y = Base::new(23);

    assert_eq!(F::from(x) * F::from(y), F::from(x * y));
}

#[test]
fn extension_frobenius_preserves_addition_and_multiplication() {
    let a = element(core::array::from_fn(|i| ((i * 3 + 2) % 31) as u8));
    let b = element(core::array::from_fn(|i| ((i * 17 + 5) % 31) as u8));

    assert_eq!((a + b).frobenius(), a.frobenius() + b.frobenius());
    assert_eq!((a * b).frobenius(), a.frobenius() * b.frobenius());
}

#[test]
fn extension_frobenius_period_is_degree() {
    let a = element(core::array::from_fn(|i| ((i * 19 + 3) % 31) as u8));

    assert_eq!(a.frobenius_pow(61), a);
}

#[test]
fn extension_generator_satisfies_configured_polynomial() {
    let xi = F::generator();
    let mut power = F::one();

    // Evaluate X^61 - (28 + 6X + 30X^2) at the canonical generator.
    for _ in 0..61 {
        power *= xi;
    }

    let rhs = F::from(Base::new(28))
        + F::from(Base::new(6)) * xi
        + F::from(Base::new(30)) * xi * xi;

    assert_eq!(power, rhs);
}

// Keep the concrete configuration in scope so this test fails to compile if
// the field alias and its reduction configuration drift apart.
const _: Option<Level1Config> = None;
