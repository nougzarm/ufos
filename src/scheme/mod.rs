use std::assert_eq;
use std::marker::PhantomData;

use crate::field::{
    extension::Fq,
    prime_field::Fp,
    traits::{ExtensionConfig, FiniteField},
};
use crate::matrix::Matrix;
use crate::random::{sample_salt, sample_secret_seed};
use crate::scheme::traits::{SignatureScheme, UfoParameters};

pub mod traits;

/// The UFOs (*Unbalanced Frobenius Oil and Vinegar*) post-quantum signature scheme.
///
/// Implements the structured UFOs signature scheme based on Frobenius automorphisms over
/// the extension field $\mathbb{F}_{p^e}$.
///
/// # Type Parameters
///
/// * `P` - Prime field characteristic ($\mathbb{F}_p$).
/// * `E` - Extension degree ($\mathbb{F}_{p^e} = \mathbb{F}_q$).
/// * `C` - Extension field configuration and exponent schedule (`ExtensionConfig + UfoParameters`).
/// * `L` - Number of Frobenius blocks per public equation.
/// * `N` - Total number of variables ($N = O + V$).
/// * `M` - Number of public equations (equal to $O$).
/// * `O` - Dimension of the secret oil subspace.
/// * `V` - Dimension of the vinegar subspace.
///
/// # Examples
///
/// ```rust
/// use ufos::params::Level1Config;
/// use ufos::scheme::UFOs;
/// use ufos::scheme::traits::SignatureScheme;
///
/// type UfosLevel1 = UFOs<31, 61, Level1Config, 5, 5, 1, 1, 4>;
/// let scheme = UfosLevel1::new();
/// let (sk, pk) = scheme.keygen();
/// ```
pub struct UFOs<
    const P: u8,
    const E: usize,
    C,
    const L: usize,
    const N: usize,
    const M: usize,
    const O: usize,
    const V: usize,
> where
    C: ExtensionConfig<E, Field = Fp<P>> + UfoParameters<L, M>,
{
    _config: PhantomData<C>,
}

impl<
        const P: u8,
        const E: usize,
        C,
        const L: usize,
        const N: usize,
        const M: usize,
        const O: usize,
        const V: usize,
    > UFOs<P, E, C, L, N, M, O, V>
where
    C: ExtensionConfig<E, Field = Fp<P>> + UfoParameters<L, M>,
{
    pub fn new() -> Self {
        Self {
            _config: PhantomData,
        }
    }

    pub fn compute_lambda_str(secret_seed: &[u8; 32]) -> (Matrix<Fq<P, E, C>>, [u8; 16]) {
        let lambda_coeffs = Fq::<P, E, C>::expand_from_seed(secret_seed, O * V);
        let lambda_str = Matrix::<Fq<P, E, C>>::new(O, V, lambda_coeffs);

        let public_seed = [0u8; 16];

        (lambda_str, public_seed)
    }

    pub fn compute_a_str(
        public_seed: &[u8; 16],
    ) -> (
        [[Matrix<Fq<P, E, C>>; L]; M],
        [[Matrix<Fq<P, E, C>>; L]; M],
        [[Matrix<Fq<P, E, C>>; L]; M],
    ) {
        let public_elements = Fq::<P, E, C>::expand_from_seed(public_seed, M * L * (N * N - O * O));

        let a_1: [[Matrix<Fq<P, E, C>>; L]; M] = core::array::from_fn(|r| {
            let a_1_r: [Matrix<Fq<P, E, C>>; L] = core::array::from_fn(|k| {
                let a_1_coeffs = &public_elements
                    [(r * L + k) * (N * N - O * O)..(r * L + k) * (N * N - O * O) + V * O];
                Matrix::<Fq<P, E, C>>::from((V, O, a_1_coeffs))
            });
            a_1_r
        });

        let a_2: [[Matrix<Fq<P, E, C>>; L]; M] = core::array::from_fn(|r| {
            let a_2_r: [Matrix<Fq<P, E, C>>; L] = core::array::from_fn(|k| {
                let a_2_coeffs = &public_elements[(r * L + k) * (N * N - O * O) + V * O
                    ..(r * L + k) * (N * N - O * O) + 2 * V * O];
                Matrix::<Fq<P, E, C>>::from((O, V, a_2_coeffs))
            });
            a_2_r
        });

        let a_3: [[Matrix<Fq<P, E, C>>; L]; M] = core::array::from_fn(|r| {
            let a_3_r: [Matrix<Fq<P, E, C>>; L] = core::array::from_fn(|k| {
                let a_3_coeffs = &public_elements
                    [(r * L + k) * (N * N - O * O) + 2 * V * O..(r * L + k + 1) * (N * N - O * O)];
                Matrix::<Fq<P, E, C>>::from((V, V, a_3_coeffs))
            });
            a_3_r
        });

        (a_1, a_2, a_3)
    }

    pub fn compute_a_zero(
        a_str: &(
            [[Matrix<Fq<P, E, C>>; L]; M],
            [[Matrix<Fq<P, E, C>>; L]; M],
            [[Matrix<Fq<P, E, C>>; L]; M],
        ),
        lambda_str: &Matrix<Fq<P, E, C>>,
    ) -> [[Matrix<Fq<P, E, C>>; L]; M] {
        let a_zero: [[Matrix<Fq<P, E, C>>; L]; M] = core::array::from_fn(|r| {
            let a_zero_r: [Matrix<Fq<P, E, C>>; L] = core::array::from_fn(|k| {
                -(&lambda_str.frobenius_pow(C::A_SCHEDULE[r][k]) * &a_str.0[r][k]
                    + &a_str.1[r][k] * &lambda_str.frobenius_pow(C::B_SCHEDULE[r][k]).transpose()
                    + &lambda_str.frobenius_pow(C::A_SCHEDULE[r][k])
                        * &a_str.2[r][k]
                        * &lambda_str.frobenius_pow(C::B_SCHEDULE[r][k]).transpose())
            });

            a_zero_r
        });

        a_zero
    }

    pub fn compute_a(
        a_str: &(
            [[Matrix<Fq<P, E, C>>; L]; M],
            [[Matrix<Fq<P, E, C>>; L]; M],
            [[Matrix<Fq<P, E, C>>; L]; M],
        ),
        a_zero: &[[Matrix<Fq<P, E, C>>; L]; M],
    ) -> [[Matrix<Fq<P, E, C>>; L]; M] {
        let a: [[Matrix<Fq<P, E, C>>; L]; M] = core::array::from_fn(|r| {
            let a_r: [Matrix<Fq<P, E, C>>; L] = core::array::from_fn(|k| {
                Matrix::<Fq<P, E, C>>::from_matrices(
                    2,
                    2,
                    &[
                        a_zero[r][k].clone(),
                        a_str.1[r][k].clone(),
                        a_str.0[r][k].clone(),
                        a_str.2[r][k].clone(),
                    ],
                )
            });
            a_r
        });
        a
    }

    pub fn la_matrix(
        lambda_str: &Matrix<Fq<P, E, C>>,
        a_2: &[[Matrix<Fq<P, E, C>>; L]; M],
        a_3: &[[Matrix<Fq<P, E, C>>; L]; M],
    ) -> [[Matrix<Fq<P, E, C>>; L]; M] {
        let la: [[Matrix<Fq<P, E, C>>; L]; M] = core::array::from_fn(|r| {
            let la_r: [Matrix<Fq<P, E, C>>; L] = core::array::from_fn(|k| {
                (&lambda_str.frobenius_pow(C::A_SCHEDULE[r][k]) * &a_3[r][k] + &a_2[r][k])
                    .transpose()
            });

            la_r
        });

        la
    }

    pub fn lb_matrix(
        lambda_str: &Matrix<Fq<P, E, C>>,
        a_1: &[[Matrix<Fq<P, E, C>>; L]; M],
        a_3: &[[Matrix<Fq<P, E, C>>; L]; M],
    ) -> [[Matrix<Fq<P, E, C>>; L]; M] {
        let lb: [[Matrix<Fq<P, E, C>>; L]; M] = core::array::from_fn(|r| {
            let lb_r: [Matrix<Fq<P, E, C>>; L] = core::array::from_fn(|k| {
                &a_3[r][k] * &lambda_str.frobenius_pow(C::B_SCHEDULE[r][k]).transpose() + &a_1[r][k]
            });

            lb_r
        });

        lb
    }
}

impl<
        const P: u8,
        const E: usize,
        C,
        const L: usize,
        const N: usize,
        const M: usize,
        const O: usize,
        const V: usize,
    > SignatureScheme for UFOs<P, E, C, L, N, M, O, V>
where
    C: ExtensionConfig<E, Field = Fp<P>> + UfoParameters<L, M>,
{
    type SecretKey = [u8; 32];
    type PublicKey = ([u8; 16], [[Matrix<Fq<P, E, C>>; L]; M]);
    type Message = Vec<u8>;
    type Signature = ([Fq<P, E, C>; N], [u8; 32]);

    fn keygen(&self) -> (Self::SecretKey, Self::PublicKey) {
        let secret_seed = sample_secret_seed();

        let (lambda_str, public_seed) = Self::compute_lambda_str(&secret_seed);

        let a_str = Self::compute_a_str(&public_seed);

        let a_zero = Self::compute_a_zero(&a_str, &lambda_str);

        (secret_seed, (public_seed, a_zero))
    }

    fn sign(&self, sk: &Self::SecretKey, m: &Self::Message) -> Self::Signature {
        let (lambda_str, public_seed) = Self::compute_lambda_str(sk);

        let a_str = Self::compute_a_str(&public_seed);

        let salt = sample_salt();

        let xi = Fq::<P, E, C>::base();

        let mut hasher = sha3::Shake256::default();
        use sha3::digest::{ExtendableOutput, Update, XofReader};
        hasher.update(m);
        hasher.update(&salt);
        let mut reader = hasher.finalize_xof();
        let mut hash_seed = [0u8; 32];
        reader.read(&mut hash_seed);

        let la_matrix = Self::la_matrix(&lambda_str, &a_str.1, &a_str.2);
        let lb_matrix = Self::lb_matrix(&lambda_str, &a_str.0, &a_str.2);

        let h_vec = Fq::<P, E, C>::expand_from_seed(&hash_seed, M);

        loop {
            let mut v_rng = rand::thread_rng();
            let mut v_coeffs = [Fq::<P, E, C>::zero(); N];
            for i in O..N {
                v_coeffs[i] = <Fq<P, E, C>>::random(&mut v_rng);
            }
            let v_str = Matrix::<Fq<P, E, C>>::from((1, V, &v_coeffs[O..N]));
            let v = Matrix::<Fq<P, E, C>>::from((1, N, &v_coeffs[0..N]));

            let mut alpha = [[[Fq::<P, E, C>::zero(); O]; L]; M];
            let mut beta = [[[Fq::<P, E, C>::zero(); O]; L]; M];
            let mut xi_pows = [[[(Fq::<P, E, C>::zero(), Fq::<P, E, C>::zero()); E]; L]; M];

            for r in 0..M {
                for k in 0..L {
                    let va_tmp = v_str.frobenius_pow(C::B_SCHEDULE[r][k]);
                    let vb_tmp = v_str.frobenius_pow(C::A_SCHEDULE[r][k]);
                    alpha[r][k] = core::array::from_fn(|j| (&va_tmp * &la_matrix[r][k])[(0, j)]);

                    beta[r][k] = core::array::from_fn(|j| (&vb_tmp * &lb_matrix[r][k])[(0, j)]);

                    xi_pows[r][k] = core::array::from_fn(|j| {
                        (
                            xi[j].frobenius_pow(C::A_SCHEDULE[r][k]),
                            xi[j].frobenius_pow(C::B_SCHEDULE[r][k]),
                        )
                    });
                }
            }

            let mut l_v_xi = [[[Fq::<P, E, C>::zero(); E]; O]; M];
            let mut pub_v = vec![Fq::<P, E, C>::zero(); M];
            let mut m_v = Matrix::<Fp<P>>::zero(O * E, M * E);
            let mut gamma = Matrix::<Fp<P>>::zero(1, M * E);

            for r in 0..M {
                for j in 0..O {
                    for u in 0..E {
                        for k in 0..L {
                            l_v_xi[r][j][u] += alpha[r][k][j] * xi_pows[r][k][u].0
                                + beta[r][k][j] * xi_pows[r][k][u].1;
                        }
                    }
                }

                for k in 0..L {
                    pub_v[r] += (&v_str.frobenius_pow(C::A_SCHEDULE[r][k])
                        * &a_str.2[r][k]
                        * &v_str.frobenius_pow(C::B_SCHEDULE[r][k]).transpose())
                        .into_scalar();
                }

                for i in 0..E {
                    gamma[(0, r * E + i)] = (h_vec[r] - pub_v[r]).coeffs[i];
                }
            }

            for r in 0..M {
                for j in 0..O {
                    for u in 0..E {
                        for i2 in 0..E {
                            m_v[(j * E + u, r * E + i2)] = l_v_xi[r][j][u].coeffs[i2];
                        }
                    }
                }
            }

            if let Some(x_fp) = m_v.solve_left(&gamma.coeffs) {
                assert_eq!(Matrix::<Fp<P>>::new(1, O * E, x_fp.clone()) * m_v, gamma);

                let x_fq_coeffs: Vec<Fq<P, E, C>> = x_fp
                    .chunks_exact(E)
                    .map(|chunk| Fq::<P, E, C>::new(chunk.try_into().unwrap()))
                    .collect();
                let x_0 = Matrix::<Fq<P, E, C>>::new(1, O, x_fq_coeffs);

                let sigma_mat =
                    x_0 * &Matrix::<Fq<P, E, C>>::from_matrices(
                        1,
                        2,
                        &[Matrix::<Fq<P, E, C>>::identity(O), lambda_str],
                    ) + v;

                let sigma: [Fq<P, E, C>; N] = core::array::from_fn(|i| sigma_mat.coeffs[i]);

                return (sigma, salt);
            }
        }
    }

    fn verify(
        &self,
        pk: &Self::PublicKey,
        m: &Self::Message,
        s: &Self::Signature,
    ) -> Result<(), ()> {
        let a_str = Self::compute_a_str(&pk.0);

        let a = Self::compute_a(&a_str, &pk.1);

        let mut hasher = sha3::Shake256::default();
        use sha3::digest::{ExtendableOutput, Update, XofReader};
        hasher.update(m);
        hasher.update(&s.1);
        let mut reader = hasher.finalize_xof();
        let mut hash_seed = [0u8; 32];
        reader.read(&mut hash_seed);

        let h_vec = Fq::<P, E, C>::expand_from_seed(&hash_seed, M);

        let mut pub_eval: Vec<Fq<P, E, C>> = vec![Fq::<P, E, C>::zero(); M];

        let sigma = Matrix::<Fq<P, E, C>>::from((1, N, &s.0[0..N]));

        for (r, item) in pub_eval.iter_mut().enumerate() {
            for k in 0..L {
                *item += (&sigma.frobenius_pow(C::A_SCHEDULE[r][k])
                    * &a[r][k]
                    * &sigma.transpose().frobenius_pow(C::B_SCHEDULE[r][k]))
                    .into_scalar();
            }

            if *item != h_vec[r] {
                return Err(());
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::extension::Fq;
    use crate::field::traits::FiniteField;
    use crate::matrix::Matrix;
    use crate::params::Level1Config;
    use crate::scheme::traits::{SignatureScheme, UfoParameters};
    use crate::UfosLevel1;

    #[test]
    fn test_oil_subspace_isotropy_constraint() {
        let scheme = UfosLevel1::new();
        let (sk, pk) = scheme.keygen();

        let (lambda_str, _) = UfosLevel1::compute_lambda_str(&sk);
        let a_str = UfosLevel1::compute_a_str(&pk.0);
        let a = UfosLevel1::compute_a(&a_str, &pk.1);

        let lambda = Matrix::<Fq<31, 61, Level1Config>>::from_matrices(
            1,
            2,
            &[Matrix::identity(1), lambda_str],
        );

        for r in 0..1 {
            for k in 0..2 {
                let a_k = Level1Config::A_SCHEDULE[r][k];
                let b_k = Level1Config::B_SCHEDULE[r][k];

                let prod = lambda.frobenius_pow(a_k)
                    * a[r][k].clone()
                    * lambda.frobenius_pow(b_k).transpose();

                assert_eq!(prod, Matrix::zero(1, 1));
            }
        }
    }

    #[test]
    fn test_level1_sign_and_verify() {
        let scheme = UfosLevel1::new();
        let (sk, pk) = scheme.keygen();
        let msg = b"Post-quantum signature for UFOs level I (p=31, e=61)".to_vec();

        let sig = scheme.sign(&sk, &msg);
        let res = scheme.verify(&pk, &msg, &sig);

        assert!(res.is_ok());
    }

    #[test]
    fn test_reject_tampered_message() {
        let scheme = UfosLevel1::new();
        let (sk, pk) = scheme.keygen();
        let msg = b"Original message".to_vec();

        let sig = scheme.sign(&sk, &msg);

        let mut wrong_msg = msg.clone();
        wrong_msg[0] ^= 0x01;

        assert!(scheme.verify(&pk, &wrong_msg, &sig).is_err());
    }

    #[test]
    fn test_reject_tampered_signature_vector() {
        let scheme = UfosLevel1::new();
        let (sk, pk) = scheme.keygen();
        let msg = b"Original message".to_vec();

        let mut sig = scheme.sign(&sk, &msg);
        // On modifie la première coordonnée du vecteur Sigma dans Fq
        sig.0[0] = sig.0[0] + Fq::<31, 61, Level1Config>::one();

        assert!(scheme.verify(&pk, &msg, &sig).is_err());
    }

    #[test]
    fn test_reject_tampered_salt() {
        let scheme = UfosLevel1::new();
        let (sk, pk) = scheme.keygen();
        let msg = b"Original message".to_vec();

        let mut sig = scheme.sign(&sk, &msg);
        sig.1[0] ^= 0xFF;

        assert!(scheme.verify(&pk, &msg, &sig).is_err());
    }

    #[test]
    fn test_reject_wrong_public_key() {
        let scheme = UfosLevel1::new();
        let (sk1, _pk1) = scheme.keygen();
        let (_sk2, pk2) = scheme.keygen();
        let msg = b"Original message".to_vec();

        // On signe avec sk1 mais on vérifie avec pk2
        let sig = scheme.sign(&sk1, &msg);

        assert!(scheme.verify(&pk2, &msg, &sig).is_err());
    }
}
