use std::ops::{Add, Sub, Mul, Div, Rem, AddAssign, SubAssign, MulAssign};
use std::fmt::Debug;
use std::hash::Hash;

use crate::my_float::MyFloatTrait;
use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::matrix_element_position::MatrixElementPosition;


pub fn compare_with_tolerance<V>(value: V, tolerance: V) -> V
    where V: MyFloatTrait + PartialOrd + From<f32>
{
    if value.my_abs() < tolerance
    {
        V::from(0f32)
    }
    else
    {
        value
    }
}


pub fn extract_unique_elements_of_rotation_matrix<T, V>(rotation_matrix: &ExtendedMatrix<T, V>)
    -> Result<Vec<V>, String>
    where T: Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> +
             Rem<Output = T> + Copy + Debug + Eq + Hash + AddAssign + SubAssign + From<u8> +
             PartialOrd + Ord + 'static,
          V: Add<Output = V> + Sub<Output = V> + Mul<Output = V> + Div<Output = V> + Copy +
             Debug + AddAssign + SubAssign + MulAssign + PartialEq + From<f32> + Into<f64> +
             'static
{
    let r11 = rotation_matrix.copy_element_value_or_zero(
        MatrixElementPosition::create(T::from(0u8), T::from(0u8)))?;
    let r12 = rotation_matrix.copy_element_value_or_zero(
        MatrixElementPosition::create(T::from(0u8), T::from(1u8)))?;
    let r13 = rotation_matrix.copy_element_value_or_zero(
        MatrixElementPosition::create(T::from(0u8), T::from(2u8)))?;
    let r21 = rotation_matrix.copy_element_value_or_zero(
        MatrixElementPosition::create(T::from(1u8), T::from(0u8)))?;
    let r22 = rotation_matrix.copy_element_value_or_zero(
        MatrixElementPosition::create(T::from(1u8), T::from(1u8)))?;
    let r23 = rotation_matrix.copy_element_value_or_zero(
        MatrixElementPosition::create(T::from(1u8), T::from(2u8)))?;
    let r31 = rotation_matrix.copy_element_value_or_zero(
        MatrixElementPosition::create(T::from(2u8), T::from(0u8)))?;
    let r32 = rotation_matrix.copy_element_value_or_zero(
        MatrixElementPosition::create(T::from(2u8), T::from(1u8)))?;
    let r33 = rotation_matrix.copy_element_value_or_zero(
        MatrixElementPosition::create(T::from(2u8), T::from(2u8)))?;
    Ok(vec![r11, r12, r13, r21, r22, r23, r31, r32, r33])
}
