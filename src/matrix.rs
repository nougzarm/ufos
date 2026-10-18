use core::ops::{Add, Mul, Sub};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Matrix<F, const R: usize, const C: usize> {
    pub coeffs: [[F; C]; R],
}

impl<F, const R: usize, const C: usize> Matrix<F, R, C> {
    fn new(coeffs: [[F; C]; R]) -> Self {
        Self { coeffs }
    }

    pub fn transpose(self) -> Matrix<F, C, R>
    where
        F: Default + Copy,
    {
        let mut transposed = [[F::default(); R]; C];

        for i in 0..R {
            for j in 0..C {
                transposed[j][i] = self.coeffs[i][j];
            }
        }

        Matrix::new(transposed)
    }
}

impl<G, F, const R: usize, const C: usize> From<[[G; C]; R]> for Matrix<F, R, C>
where
    G: Copy,
    F: From<G>,
{
    fn from(value: [[G; C]; R]) -> Self {
        Self::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| F::from(value[i][j]))
        }))
    }
}

impl<F, const R: usize, const C: usize> Add for Matrix<F, R, C>
where
    F: Add<Output = F> + Copy,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| self.coeffs[i][j] + rhs.coeffs[i][j])
        }))
    }
}

impl<F, const R: usize, const C: usize> Sub for Matrix<F, R, C>
where
    F: Sub<Output = F> + Copy,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| self.coeffs[i][j] - rhs.coeffs[i][j])
        }))
    }
}

impl<F, const R1: usize, const K: usize, const C2: usize> Mul<Matrix<F, K, C2>> for Matrix<F, R1, K>
where
    F: Add<Output = F> + Mul<Output = F> + Default + Copy,
{
    type Output = Matrix<F, R1, C2>;

    fn mul(self, rhs: Matrix<F, K, C2>) -> Self::Output {
        let mut result = [[F::default(); C2]; R1];

        for i in 0..R1 {
            for j in 0..C2 {
                for k in 0..K {
                    result[i][j] = result[i][j] + self.coeffs[i][k] * rhs.coeffs[k][j];
                }
            }
        }

        Matrix::new(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::prime_field::Fp;

    type F31 = Fp<31>;

    #[test]
    fn test_matrix_add_sub() {
        let a = Matrix::<F31, 2, 2>::from([[1, 2], [3, 4]]);

        let b = Matrix::<F31, 2, 2>::from([[10, 20], [30, 9]]);

        let c = a + b;
        assert_eq!(c, Matrix::<F31, 2, 2>::from([[11, 22], [2, 13],]));

        let d = b - a;
        assert_eq!(d, Matrix::<F31, 2, 2>::from([[9, 18], [27, 5],]));
    }

    #[test]
    fn test_matrix_rect_mul() {
        let a = Matrix::<F31, 2, 3>::from([[1, 2, 3], [4, 5, 6]]);

        let b = Matrix::<F31, 3, 2>::from([[7, 8], [9, 1], [2, 3]]);

        let c = a * b;
        assert_eq!(c, Matrix::<F31, 2, 2>::from([[0, 19], [85, 55],]));
    }
}
