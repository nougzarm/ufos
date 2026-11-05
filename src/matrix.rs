//! # Linear Algebra and Matrix Operations
//!
//! This module provides the [`Matrix`] struct for performing linear algebra
//! operations over any field implementing the [`FiniteField`] trait.
//!
//! Key features include:
//! - Block-matrix assembly via [`Matrix::from_matrices`].
//! - Linear system solving over base fields via [`Matrix::solve_left`].
//! - Frobenius endomorphism application via [`Matrix::frobenius_pow`].

use core::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::field::{
    extension::Fq,
    prime_field::Fp,
    traits::{ExtensionConfig, FiniteField},
};

/// A structure representing a matrix with a generic entries.
///
/// Elements are stored sequentially in a single contiguous [`Vec`] indexed as
/// `(i, j) \mapsto i \times \mathrm{cols} + j$.
///
/// # Type Parameters
///
/// * `F` - Scalar structure of the matrix entries.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<F> {
    pub rows: usize,
    pub cols: usize,
    pub coeffs: Vec<F>,
}

impl<F> Index<(usize, usize)> for Matrix<F> {
    type Output = F;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        debug_assert!(i < self.rows && j < self.cols);
        &self.coeffs[i * self.cols + j]
    }
}

impl<F> IndexMut<(usize, usize)> for Matrix<F> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        debug_assert!(i < self.rows && j < self.cols);
        &mut self.coeffs[i * self.cols + j]
    }
}

impl<G, F, const R: usize, const C: usize> From<[[G; C]; R]> for Matrix<F>
where
    G: Copy,
    F: From<G>,
{
    fn from(value: [[G; C]; R]) -> Self {
        Self {
            rows: R,
            cols: C,
            coeffs: value
                .as_flattened()
                .to_vec()
                .into_iter()
                .map(Into::into)
                .collect(),
        }
    }
}

impl<F> From<(usize, usize, &[F])> for Matrix<F>
where
    F: Copy,
{
    fn from(value: (usize, usize, &[F])) -> Self {
        Self::new(value.0, value.1, value.2.to_vec())
    }
}

impl<F> Matrix<F>
where
    F: Copy,
{
    /// Creates a new matrix from row and column dimensions and a flattened coefficient buffer.
    ///
    /// # Arguments
    ///
    /// * `rows` - Number of rows.
    /// * `cols` - Number of columns.
    /// * `coeffs` - Row-major vector containing exactly `rows * cols` entries.
    pub fn new(rows: usize, cols: usize, coeffs: Vec<F>) -> Self {
        Self { rows, cols, coeffs }
    }

    /// Assembles a global matrix from a two-dimensional grid of sub-matrix blocks.
    ///
    /// Given a grid of dimensions $B_r \times B_c$, this function concatenates the constituent
    /// blocks row-by-row into a single dense matrix.
    ///
    /// # Arguments
    ///
    /// * `blocs_rows` - Number of block rows ($B_r$).
    /// * `blocs_cols` - Number of block columns ($B_c$).
    /// * `matrices` - Flat slice of blocks in row-major order, with length $B_r \times B_c$.
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - `blocs_rows == 0` or `blocs_cols == 0`.
    /// - `matrices.len() != blocs_rows * blocs_cols`.
    /// - Sub-matrices in the same row do not share identical row heights.
    /// - Sub-matrices in the same column do not share identical column widths.
    pub fn from_matrices(blocs_rows: usize, blocs_cols: usize, matrices: &[Matrix<F>]) -> Self {
        assert_eq!(blocs_rows * blocs_cols, matrices.len());
        assert!(blocs_rows > 0 && blocs_cols > 0);

        let mut row_heights = vec![0; blocs_rows];
        let mut col_widths = vec![0; blocs_cols];

        for r in 0..blocs_rows {
            for c in 0..blocs_cols {
                let idx = r * blocs_cols + c;
                let mat = &matrices[idx];

                if r == 0 {
                    col_widths[c] = mat.cols;
                } else {
                    assert_eq!(col_widths[c], mat.cols);
                }

                if c == 0 {
                    row_heights[r] = mat.rows;
                } else {
                    assert_eq!(row_heights[r], mat.rows);
                }
            }
        }

        let total_rows: usize = row_heights.iter().sum();
        let total_cols: usize = col_widths.iter().sum();

        let mut res_coeffs = Vec::with_capacity(total_rows * total_cols);
        for r in 0..blocs_rows {
            let h = row_heights[r];
            for local_i in 0..h {
                for c in 0..blocs_cols {
                    let idx = r * blocs_cols + c;
                    let mat = &matrices[idx];
                    let w = col_widths[c];

                    let src_start = local_i * w;
                    let src_end = src_start + w;
                    res_coeffs.extend_from_slice(&mat.coeffs[src_start..src_end]);
                }
            }
        }

        Self {
            rows: total_rows,
            cols: total_cols,
            coeffs: res_coeffs,
        }
    }

    /// Computes the transpose of a matrix.
    ///
    /// For $M \in \mathbb{F}^{R \times C}$, returns ${}^t M \in \mathbb{F}^{C \times R}$
    /// satisfying $({}^t M)_{j, i} = M_{i, j}$.
    pub fn transpose(&self) -> Self {
        let mut res = Self {
            rows: self.cols,
            cols: self.rows,
            coeffs: self.coeffs.clone(),
        };
        for i in 0..self.rows {
            for j in 0..self.cols {
                res[(j, i)] = self[(i, j)];
            }
        }
        res
    }

    /// Extracts the single scalar entry of a $1 \times 1$ matrix.
    ///
    /// # Panics
    ///
    /// Panics if the matrix dimension is not $1 \times 1$.
    pub fn into_scalar(&self) -> F {
        assert_eq!((self.rows, self.cols), (1, 1));
        self[(0, 0)]
    }
}

impl<F> Matrix<F>
where
    F: FiniteField,
{
    /// Instantiates an all-zero matrix $0_{rows \times cols} \in \mathbb{F}^{rows \times cols}$.
    ///
    /// # Arguments
    ///
    /// * `rows` - Number of rows.
    /// * `cols` - Number of columns.
    pub fn zero(rows: usize, cols: usize) -> Self {
        let coeffs = vec![F::zero(); rows * cols];
        Self::new(rows, cols, coeffs)
    }

    /// Generates the $n \times n$ identity matrix $1_n$.
    ///
    /// # Arguments
    ///
    /// * `dimension` - Size $n$ of the square matrix.
    pub fn identity(dimension: usize) -> Self {
        let mut res = Self::zero(dimension, dimension);
        for i in 0..dimension {
            res[(i, i)] = F::one();
        }

        res
    }

    /// Solves the left linear system $X \cdot M = Y$ over a finite field $\mathbb{F}$ using Gauss-Jordan elimination.
    ///
    /// Finds the row vector $X \in \mathbb{F}^{1 \times n}$ satisfying the matrix equation
    /// for a square matrix $A \in \mathbb{F}^{n \times n}$ and right-hand side $Y \in \mathbb{F}^{1 \times n}$.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The right-hand side vector $Y$, with length equal to `self.rows`.
    ///
    /// # Returns
    ///
    /// Returns `Some(X)` if the system has a unique solution, or `None` if the
    /// matrix is singular ($\det(M) = 0$).
    ///
    /// # Panics
    ///
    /// Panics if the matrix is not square (`self.rows != self.cols`) or if the
    /// length of `rhs` does not match the matrix dimension.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use ufos::matrix::Matrix;
    /// use ufos::field::prime_field::Fp;
    ///
    /// let id = Matrix::<Fp<31>>::identity(3);
    /// let rhs = vec![Fp::new(1), Fp::new(2), Fp::new(3)];
    /// let sol = id.solve_left(&rhs).unwrap();
    /// assert_eq!(sol, rhs);
    /// ```
    pub fn solve_left(&self, rhs: &[F]) -> Option<Vec<F>> {
        let n: usize = self.rows;
        assert_eq!(n, self.cols);
        assert_eq!(n, rhs.len());

        let mut a = self.transpose();
        let mut b = rhs.to_vec();

        for i in 0..n {
            let mut pivot = i;
            while pivot < n && a[(pivot, i)] == F::zero() {
                pivot += 1;
            }
            if pivot == n {
                return None;
            }

            if pivot != i {
                for col in 0..n {
                    let tmp = a[(i, col)];
                    a[(i, col)] = a[(pivot, col)];
                    a[(pivot, col)] = tmp;
                }
                b.swap(i, pivot);
            }

            let inv_pivot = a[(i, i)].inv();
            for col in i..n {
                a[(i, col)] *= inv_pivot;
            }
            b[i] *= inv_pivot;

            for row in 0..n {
                if row != i {
                    let factor = a[(row, i)];
                    if factor != F::zero() {
                        for col in i..n {
                            a[(row, col)] = a[(row, col)] - factor * a[(i, col)];
                        }
                        b[row] = b[row] - factor * b[i];
                    }
                }
            }
        }

        Some(b)
    }
}

impl<const P: u8, const E: usize, C> Matrix<Fq<P, E, C>>
where
    C: ExtensionConfig<E, Field = Fp<P>>,
{
    /// Applies the $a$-th Frobenius automorphism component-wise to all matrix entries.
    ///
    /// Computes $M^{[a]} = (m_{i, j}^{p^a})$ for each entry $m_{i, j} \in \mathbb{F}_{p^e}$,
    /// where exponents $a$ are interpreted modulo $e$.
    ///
    /// # Arguments
    ///
    /// * `a` - Exponent power of the Frobenius automorphism $\phi^a : x \mapsto x^{p^a}$.
    pub fn frobenius_pow(&self, a: usize) -> Self {
        let mut res = self.clone();

        for c in res.coeffs.iter_mut() {
            *c = c.frobenius_pow(a);
        }

        res
    }
}

impl<F> Add<&Matrix<F>> for &Matrix<F>
where
    F: Add<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn add(self, rhs: &Matrix<F>) -> Self::Output {
        assert_eq!((self.rows, self.cols), (rhs.rows, rhs.cols));
        let coeffs = self
            .coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .map(|(&a, &b)| a + b)
            .collect();
        Matrix {
            rows: self.rows,
            cols: self.cols,
            coeffs,
        }
    }
}

impl<F> Add<Matrix<F>> for Matrix<F>
where
    F: Add<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn add(self, rhs: Matrix<F>) -> Self::Output {
        &self + &rhs
    }
}

impl<F> Add<&Matrix<F>> for Matrix<F>
where
    F: Add<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn add(self, rhs: &Matrix<F>) -> Self::Output {
        &self + rhs
    }
}

impl<F> Add<Matrix<F>> for &Matrix<F>
where
    F: Add<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn add(self, rhs: Matrix<F>) -> Self::Output {
        self + &rhs
    }
}

impl<F> AddAssign<&Matrix<F>> for Matrix<F>
where
    F: Add<Output = F> + Copy,
{
    fn add_assign(&mut self, rhs: &Matrix<F>) {
        assert_eq!((self.rows, self.cols), (rhs.rows, rhs.cols),);
        for (a, &b) in self.coeffs.iter_mut().zip(rhs.coeffs.iter()) {
            *a = *a + b;
        }
    }
}

impl<F> AddAssign<Matrix<F>> for Matrix<F>
where
    F: Add<Output = F> + Copy,
{
    fn add_assign(&mut self, rhs: Matrix<F>) {
        *self += &rhs;
    }
}

impl<F> Sub<&Matrix<F>> for &Matrix<F>
where
    F: Sub<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn sub(self, rhs: &Matrix<F>) -> Self::Output {
        assert_eq!((self.rows, self.cols), (rhs.rows, rhs.cols));
        let coeffs = self
            .coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .map(|(&a, &b)| a - b)
            .collect();
        Matrix {
            rows: self.rows,
            cols: self.cols,
            coeffs,
        }
    }
}

impl<F> Sub<Matrix<F>> for Matrix<F>
where
    F: Sub<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn sub(self, rhs: Matrix<F>) -> Self::Output {
        &self - &rhs
    }
}

impl<F> Sub<&Matrix<F>> for Matrix<F>
where
    F: Sub<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn sub(self, rhs: &Matrix<F>) -> Self::Output {
        &self - rhs
    }
}

impl<F> Sub<Matrix<F>> for &Matrix<F>
where
    F: Sub<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn sub(self, rhs: Matrix<F>) -> Self::Output {
        self - &rhs
    }
}

impl<F> SubAssign<&Matrix<F>> for Matrix<F>
where
    F: Sub<Output = F> + Copy,
{
    fn sub_assign(&mut self, rhs: &Matrix<F>) {
        assert_eq!((self.rows, self.cols), (rhs.rows, rhs.cols),);
        for (a, &b) in self.coeffs.iter_mut().zip(rhs.coeffs.iter()) {
            *a = *a - b;
        }
    }
}

impl<F> SubAssign<Matrix<F>> for Matrix<F>
where
    F: Sub<Output = F> + Copy,
{
    fn sub_assign(&mut self, rhs: Matrix<F>) {
        *self -= &rhs;
    }
}

impl<F> Neg for &Matrix<F>
where
    F: Neg<Output = F> + Copy,
{
    type Output = Matrix<F>;

    fn neg(self) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|&x| -x).collect();
        Matrix {
            rows: self.rows,
            cols: self.cols,
            coeffs,
        }
    }
}

impl<F> Neg for Matrix<F>
where
    F: Neg<Output = F> + Copy,
{
    type Output = Matrix<F>;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<F> Mul<&Matrix<F>> for &Matrix<F>
where
    F: Add<Output = F> + Mul<Output = F> + Default + Copy,
{
    type Output = Matrix<F>;

    /// Multiplies two matrices $A \cdot B$ over $\mathbb{F}$.
    ///
    /// # Panics
    ///
    /// Panics if the inner dimensions mismatch (`self.cols != rhs.rows`).
    fn mul(self, rhs: &Matrix<F>) -> Self::Output {
        assert_eq!(self.cols, rhs.rows);
        let mut coeffs = vec![F::default(); self.rows * rhs.cols];

        for i in 0..self.rows {
            let res_row_offset = i * rhs.cols;
            for k in 0..self.cols {
                let aik = self[(i, k)];
                let rhs_row_offset = k * rhs.cols;
                for j in 0..rhs.cols {
                    coeffs[res_row_offset + j] =
                        coeffs[res_row_offset + j] + aik * rhs.coeffs[rhs_row_offset + j];
                }
            }
        }
        Matrix {
            rows: self.rows,
            cols: rhs.cols,
            coeffs,
        }
    }
}

impl<F> Mul<Matrix<F>> for Matrix<F>
where
    F: Add<Output = F> + Mul<Output = F> + Default + Copy,
{
    type Output = Matrix<F>;

    fn mul(self, rhs: Matrix<F>) -> Self::Output {
        &self * &rhs
    }
}

impl<F> Mul<&Matrix<F>> for Matrix<F>
where
    F: Add<Output = F> + Mul<Output = F> + Default + Copy,
{
    type Output = Matrix<F>;

    fn mul(self, rhs: &Matrix<F>) -> Self::Output {
        &self * rhs
    }
}

impl<F> Mul<Matrix<F>> for &Matrix<F>
where
    F: Add<Output = F> + Mul<Output = F> + Default + Copy,
{
    type Output = Matrix<F>;

    fn mul(self, rhs: Matrix<F>) -> Self::Output {
        self * &rhs
    }
}

impl<F> Mul<F> for &Matrix<F>
where
    F: Mul<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn mul(self, scalar: F) -> Self::Output {
        let coeffs: Vec<F> = self
            .coeffs
            .clone()
            .into_iter()
            .map(|x| x * scalar)
            .collect();
        Matrix {
            rows: self.rows,
            cols: self.cols,
            coeffs,
        }
    }
}

impl<F> Mul<F> for Matrix<F>
where
    F: Mul<Output = F> + Copy,
{
    type Output = Matrix<F>;
    fn mul(self, scalar: F) -> Self::Output {
        &self * scalar
    }
}

impl<F> MulAssign<F> for Matrix<F>
where
    F: Mul<Output = F> + Copy,
{
    fn mul_assign(&mut self, scalar: F) {
        for x in self.coeffs.iter_mut() {
            *x = *x * scalar;
        }
    }
}

#[cfg(test)]
mod tests {
    use std::assert_eq;

    use super::*;
    use crate::field::prime_field::Fp;

    const P: u8 = 31;
    const N: usize = 5;
    type F31 = Fp<P>;

    #[test]
    fn test_matrix_add_sub() {
        let a = Matrix::<F31>::from([[1, 2], [3, 4]]);

        let b = Matrix::<F31>::from([[10, 20], [30, 9]]);

        let c = a.clone() + b.clone();
        assert_eq!(c, Matrix::<F31>::from([[11, 22], [2, 13]]));

        let d = b - a;
        assert_eq!(d, Matrix::<F31>::from([[9, 18], [27, 5]]));
    }

    #[test]
    fn test_matrix_rect_mul() {
        let a = Matrix::<F31>::from([[1, 2, 3], [4, 5, 6]]);

        let b = Matrix::<F31>::from([[7, 8], [9, 1], [2, 3]]);

        let c = a * b;
        assert_eq!(c, Matrix::<F31>::from([[0, 19], [85, 55]]));
    }

    #[test]
    fn test_matrix_from_matrices_blocks() {
        let a0 = Matrix::<F31>::from([[5]]);
        let a2 = Matrix::<F31>::from([[1, 2, 3, 4]]);
        let a1 = Matrix::<F31>::from([[10], [11], [12], [13]]);
        let a3 = Matrix::<F31>::from([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);

        let global_mat = Matrix::from_matrices(2, 2, &[a0, a2, a1, a3]);

        assert_eq!(global_mat.rows, 5);
        assert_eq!(global_mat.cols, 5);
        assert_eq!(global_mat[(0, 0)], F31::new(5));
        assert_eq!(global_mat[(0, 4)], F31::new(4));
        assert_eq!(global_mat[(4, 0)], F31::new(13));
        assert_eq!(global_mat[(4, 4)], F31::new(1));
    }

    #[test]
    fn test_pivot_gauss() {
        let x = F31::expand_from_seed(&[0u8; 16], N);
        let identity = Matrix::<F31>::identity(N);

        let x_res = identity.solve_left(x.as_slice()).unwrap();

        assert_eq!(x, x_res);
    }
}
