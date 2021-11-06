use std::ops::{Add, Sub, Mul, Div, Rem, AddAssign, SubAssign, MulAssign};
use std::fmt::Debug;
use std::hash::Hash;

use crate::my_float::MyFloatTrait;
use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::matrix_element_value_extractor;


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
             PartialOrd + 'static
{
    let r11 = matrix_element_value_extractor(T::from(0u8), T::from(0u8), &rotation_matrix)?;
    let r12 = matrix_element_value_extractor(T::from(0u8), T::from(1u8), &rotation_matrix)?;
    let r13 = matrix_element_value_extractor(T::from(0u8), T::from(2u8), &rotation_matrix)?;
    let r21 = matrix_element_value_extractor(T::from(1u8), T::from(0u8), &rotation_matrix)?;
    let r22 = matrix_element_value_extractor(T::from(1u8), T::from(1u8), &rotation_matrix)?;
    let r23 = matrix_element_value_extractor(T::from(1u8), T::from(2u8), &rotation_matrix)?;
    let r31 = matrix_element_value_extractor(T::from(2u8), T::from(0u8), &rotation_matrix)?;
    let r32 = matrix_element_value_extractor(T::from(2u8), T::from(1u8), &rotation_matrix)?;
    let r33 = matrix_element_value_extractor(T::from(2u8), T::from(2u8), &rotation_matrix)?;
    Ok(vec![r11, r12, r13, r21, r22, r23, r31, r32, r33])
}
