use std::ops::{Div, Rem, Mul, Add, AddAssign, Sub, SubAssign, MulAssign};
use std::fmt::Debug;
use std::hash::Hash;
use std::collections::HashMap;

use extended_matrix::matrix_element_position::MatrixElementPosition;
use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::conversion_uint_into_usize;


use crate::fem::global_analysis::fe_stiffness::
{
    stiffness_types_number, StiffnessGroup, StiffnessType, StiffnessGroupKey
};
use crate::fem::global_analysis::fe_dof_parameter_data::global_dof;

use crate::fem::separated_matrix::SeparatedMatrix;


pub fn add_new_stiffness_sub_groups<'a, T>(
    stiffness_groups: &mut HashMap<StiffnessGroupKey<T>, Vec<MatrixElementPosition<T>>>,
    global_group_position: T, global_group_columns_number: T, global_number_1: T,
    global_number_2: T) -> Result<(), &'a str>
    where T: Copy + Debug + Div<Output = T> + Rem<Output = T> + Mul<Output = T> + Add<Output = T> +
             PartialOrd + AddAssign + From<u8> + Eq + Hash + SubAssign
{
    let row = global_group_position / global_group_columns_number;
    let column = global_group_position % global_group_columns_number;

    let mut k = T::from(0u8);
    while k < stiffness_types_number()
    {
        let start_row = row * global_dof::<T>();

        let row_shift_init = k / T::from(2u8) * (global_dof::<T>() / T::from(2u8));

        let row_shift_final = k / T::from(2u8) * (global_dof::<T>() / T::from(2u8)) +
            (global_dof::<T>() / T::from(2u8));

        let start_column = column * global_dof::<T>();

        let column_shift_init = k % T::from(2u8) * (global_dof::<T>() / T::from(2u8));

        let column_shift_final = k % T::from(2u8) * (global_dof::<T>() / T::from(2u8)) +
            (global_dof::<T>() / T::from(2u8));

        let mut element_positions = Vec::new();
        let mut current_row = start_row + row_shift_init;
        while current_row < start_row + row_shift_final
        {
            let mut current_column  = start_column + column_shift_init;
            while current_column < start_column + column_shift_final
            {
                element_positions.push(
                    MatrixElementPosition::create(current_row, current_column));
                current_column += T::from(1u8);
            }
            current_row += T::from(1u8);
        }
        let converted_index = conversion_uint_into_usize(k);
        let stiffness_type = StiffnessType::iterator()
            .nth(converted_index)
            .ok_or("FEModel: Stiffness type could not be defined")?;
        let stiffness_group_key = StiffnessGroupKey { stiffness_type: *stiffness_type,
            number_1: global_number_1, number_2: global_number_2 };
        stiffness_groups.insert(stiffness_group_key, element_positions);
        k += T::from(1u8);
    }
    Ok(())
}


pub fn separate<T, V>(matrix: ExtendedMatrix<T, V>, positions: Vec<MatrixElementPosition<T>>,
    tolerance: V) -> Result<SeparatedMatrix<T, V>, String>
    where T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Copy + Debug + Eq + Hash + SubAssign + PartialOrd + AddAssign +
             From<u8> + Ord + 'static,
          V: Add<Output = V> + Mul<Output = V> + Sub<Output = V> + Div<Output = V> + Copy + Debug +
             PartialEq + AddAssign + MulAssign + SubAssign + Into<f64> + From<f32> + PartialOrd +
             'static
{
    let shape = matrix.copy_shape();

    let mut converted_positions_length = T::from(0u8);
    (0..positions.len()).for_each(|_| converted_positions_length += T::from(1u8));

    let k_aa_rows_number = shape.0 - converted_positions_length;

    let k_aa_columns_number = shape.1 - converted_positions_length;

    let mut k_aa_elements = Vec::new();

    let mut i = T::from(0u8);
    while i < shape.0
    {
        let mut j = T::from(0u8);
        while j < shape.1
        {
            if positions.iter().position(|p| *p.ref_row() == i).is_none() &&
                positions.iter().position(|p| *p.ref_column() == j).is_none()
            {
                let value = matrix.copy_element_value_or_zero(
                    MatrixElementPosition::create(i, j))?;
                k_aa_elements.push(value);
            }
            j += T::from(1u8);
        }
        i += T::from(1u8);
    }

    let k_aa_matrix = ExtendedMatrix::create(k_aa_rows_number,
        k_aa_columns_number, k_aa_elements, tolerance)?;

    let k_ab_rows_number = shape.0 - converted_positions_length;

    let k_ab_columns_number = converted_positions_length;

    let mut k_ab_elements = Vec::new();

    let mut i = T::from(0u8);
    while i < shape.0
    {
        if positions.iter().position(|p| *p.ref_row() == i).is_none()
        {
            for j in 0..positions.len()
            {
                let row = i;
                let column = positions[j].ref_column();
                if *column > shape.1
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kab \
                        could not be composed!".to_string());
                }
                let value = matrix.copy_element_value_or_zero(
                    MatrixElementPosition::create(row, *column))?;
                k_ab_elements.push(value);
            }
        }
        i += T::from(1u8);
    }

    let k_ab_matrix = ExtendedMatrix::create(k_ab_rows_number,
        k_ab_columns_number, k_ab_elements, tolerance)?;

    let k_ba_rows_number = converted_positions_length;

    let k_ba_columns_number = shape.1 - converted_positions_length;

    let mut k_ba_elements = Vec::new();


    for i in 0..positions.len()
    {
        let mut j = T::from(0u8);
        while j < shape.1
        {
            if positions.iter().position(|p| *p.ref_column() == j).is_none()
            {
                let row = positions[i].ref_row();
                let column = j;
                if *row > shape.0
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kba \
                        could not be composed!".to_string());
                }
                let value = matrix.copy_element_value_or_zero(
                    MatrixElementPosition::create(*row, column))?;
                k_ba_elements.push(value);
            }
            j += T::from(1u8);
        }
    }

    let k_ba_matrix = ExtendedMatrix::create(k_ba_rows_number,
        k_ba_columns_number, k_ba_elements, tolerance)?;

    let k_bb_rows_number = converted_positions_length;

    let k_bb_columns_number = converted_positions_length;

    let mut k_bb_elements = Vec::new();

    for i in 0..positions.len()
    {
        for j in 0..positions.len()
        {
            let row = positions[i].ref_row();
            let column = positions[j].ref_column();
            if *row > shape.0 || *column > shape.1
            {
                return Err("Extended matrix: Matrix could not be separated! Matrix Kbb could \
                    not be composed!".to_string());
            }
            let value = matrix.copy_element_value_or_zero(
                MatrixElementPosition::create(*row, *column))?;
            k_bb_elements.push(value);
        }
    }
    let k_bb_matrix = ExtendedMatrix::create(k_bb_rows_number,
        k_bb_columns_number, k_bb_elements, tolerance)?;

    let separated_matrix = SeparatedMatrix::create(k_aa_matrix,
        k_ab_matrix, k_ba_matrix, k_bb_matrix);

    Ok(separated_matrix)
}
