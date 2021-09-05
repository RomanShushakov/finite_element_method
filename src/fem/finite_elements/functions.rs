use std::ops::{Add, Sub, Mul, Div, Rem, AddAssign, SubAssign, MulAssign};
use std::fmt::Debug;
use std::hash::Hash;

use crate::my_float::MyFloatTrait;
use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::copy_element_value;


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
    -> Vec<V>
    where T: Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> +
             Rem<Output = T> + Copy + Debug + Eq + Hash + AddAssign + SubAssign + From<u8> +
             PartialOrd + 'static,
          V: Add<Output = V> + Sub<Output = V> + Mul<Output = V> + Div<Output = V> + Copy +
             Debug + AddAssign + SubAssign + MulAssign + PartialEq + From<f32> + Into<f64> +
             'static
{
    let all_elements_values_of_rotation_matrix =
        rotation_matrix.copy_all_elements_values();
    let r11 = copy_element_value(T::from(0u8), T::from(0u8),
        &all_elements_values_of_rotation_matrix);
    let r12 = copy_element_value(T::from(0u8), T::from(1u8),
        &all_elements_values_of_rotation_matrix);
    let r13 = copy_element_value(T::from(0u8), T::from(2u8),
        &all_elements_values_of_rotation_matrix);
    let r21 = copy_element_value(T::from(1u8), T::from(0u8),
        &all_elements_values_of_rotation_matrix);
    let r22 = copy_element_value(T::from(1u8), T::from(1u8),
        &all_elements_values_of_rotation_matrix);
    let r23 = copy_element_value(T::from(1u8), T::from(2u8),
        &all_elements_values_of_rotation_matrix);
    let r31 = copy_element_value(T::from(2u8), T::from(0u8),
        &all_elements_values_of_rotation_matrix);
    let r32 = copy_element_value(T::from(2u8), T::from(1u8),
        &all_elements_values_of_rotation_matrix);
    let r33 = copy_element_value(T::from(2u8), T::from(2u8),
        &all_elements_values_of_rotation_matrix);
    vec![r11, r12, r13, r21, r22, r23, r31, r32, r33]
}
