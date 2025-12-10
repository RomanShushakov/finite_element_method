#![allow(unused_imports)]

use extended_matrix::{BasicOperationsTrait, Position, SquareMatrix, Vector};

use crate::fem::{
    linalg::{axpy, dot, norm2, scale},
    preconditioners::apply_jacobi_preconditioner,
};

const ABS_TOL: f64 = 1e-12;

fn vecf(values: &[f64]) -> Vector<f64> {
    Vector::create(values)
}

#[test]
fn test_dot() {
    let x = vecf(&[1.0, 2.0, 3.0]);
    let y = vecf(&[4.0, 5.0, 6.0]);
    let d = dot(&x, &y).unwrap();
    assert!((d - 32.0).abs() < 1e-12);
}

#[test]
fn test_norm() {
    let x = vecf(&[3.0, 4.0]);
    let n = norm2(&x).unwrap();
    assert!((n - 5.0).abs() < 1e-12);
}

#[test]
fn test_axpy() {
    let x = vecf(&[3.0, 4.0]);
    let mut y = vecf(&[1.0, 1.0]);

    axpy(&mut y, 2.0_f64, &x).unwrap();

    assert!((y.get_element_value(&Position(0, 0)).unwrap() - 7.0).abs() < ABS_TOL);
    assert!((y.get_element_value(&Position(1, 0)).unwrap() - 9.0).abs() < ABS_TOL);
}

#[test]
fn test_scale() {
    let mut x = vecf(&[2.0, -1.0, 0.5]);
    scale(&mut x, 4.0_f64).unwrap();

    assert!((x.get_element_value(&Position(0, 0)).unwrap() - 8.0).abs() < ABS_TOL);
    assert!((x.get_element_value(&Position(1, 0)).unwrap() + 4.0).abs() < ABS_TOL);
    assert!((x.get_element_value(&Position(2, 0)).unwrap() - 2.0).abs() < ABS_TOL);
}

fn mat2x2(a11: f64, a12: f64, a21: f64, a22: f64) -> SquareMatrix<f64> {
    let mut m = SquareMatrix::create(2, &[a11, a12, a21, a22]);
    *m.get_mut_element_value(&Position(0, 0)).unwrap() = a11;
    *m.get_mut_element_value(&Position(0, 1)).unwrap() = a12;
    *m.get_mut_element_value(&Position(1, 0)).unwrap() = a21;
    *m.get_mut_element_value(&Position(1, 1)).unwrap() = a22;
    m
}

#[test]
fn test_jacobi_preconditioner_simple() {
    let a = mat2x2(2.0, 1.0, 1.0, 3.0);
    let r = vecf(&[4.0, 5.0]);

    let z = apply_jacobi_preconditioner(&a, &r).unwrap();

    let z0 = z.get_element_value(&Position(0, 0)).unwrap();
    let z1 = z.get_element_value(&Position(1, 0)).unwrap();

    assert!((z0 - 2.0).abs() < ABS_TOL);
    assert!((z1 - (5.0 / 3.0)).abs() < ABS_TOL);
}
