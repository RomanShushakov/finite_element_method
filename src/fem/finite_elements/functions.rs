use extended_matrix::{Matrix, FloatTrait, BasicOperationsTrait, Position};


pub fn compare_with_tolerance<V>(value: V, abs_tol: V) -> V
    where V: FloatTrait<Output = V>
{
    if value.my_abs() < abs_tol
    {
        V::from(0f32)
    }
    else
    {
        value
    }
}


pub fn extract_unique_elements_of_rotation_matrix<V>(rotation_matrix: &Matrix<V>) -> Result<Vec<V>, String>
    where V: FloatTrait<Output = V>
{
    let r11 = rotation_matrix.get_element_value(&Position(0, 0))?;
    let r12 = rotation_matrix.get_element_value(&Position(0, 1))?;
    let r13 = rotation_matrix.get_element_value(&Position(0, 2))?;
    let r21 = rotation_matrix.get_element_value(&Position(1, 0))?;
    let r22 = rotation_matrix.get_element_value(&Position(1, 1))?;
    let r23 = rotation_matrix.get_element_value(&Position(1, 2))?;
    let r31 = rotation_matrix.get_element_value(&Position(2, 0))?;
    let r32 = rotation_matrix.get_element_value(&Position(2, 1))?;
    let r33 = rotation_matrix.get_element_value(&Position(2, 2))?;
    Ok(vec![r11, r12, r13, r21, r22, r23, r31, r32, r33])
}
