use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::matrix_element_value_extractor;
use extended_matrix::traits::{UIntTrait, FloatTrait};


pub fn compare_with_tolerance<V>(value: V, tolerance: V) -> V
    where V: FloatTrait<Output = V, Other = V>
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
    where T: UIntTrait<Output = T>,
          V: FloatTrait<Output = V, Other = V>
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
