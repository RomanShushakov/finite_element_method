use std::ops::{Div, Rem, Mul, Add, AddAssign, Sub, SubAssign, MulAssign};
use std::fmt::Debug;
use std::hash::Hash;

use extended_matrix::basic_matrix::basic_matrix::MatrixElementPosition;
use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::{extract_element_value, conversion_uint_into_usize};


use crate::fem::global_analysis::fe_stiffness::
{
    stiffness_types_number, StiffnessGroup, StiffnessType
};
use crate::fem::global_analysis::fe_dof_parameter_data::global_dof;

use crate::fem::fe_model::SeparatedMatrix;

use crate::TOLERANCE;


pub fn compose_stiffness_sub_groups<'a, T>(global_group_position: T,
    global_group_columns_number: T, global_number_1: T, global_number_2: T)
    -> Result<Vec<StiffnessGroup<T>>, &'a str>
    where T: Copy + Debug + Div<Output = T> + Rem<Output = T> + Mul<Output = T> + Add<Output = T> +
             PartialOrd + AddAssign + From<u8>
{
    let mut stiffness_sub_groups = Vec::new();
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
        let stiffness_sub_group = StiffnessGroup { stiffness_type: *stiffness_type,
            number_1: global_number_1,
            number_2: global_number_2,
            positions: element_positions,
        };
        stiffness_sub_groups.push(stiffness_sub_group);
        k += T::from(1u8);
    }
    Ok(stiffness_sub_groups)
}


pub fn separate<'a, T, V>(matrix: ExtendedMatrix<T, V>, positions: Vec<MatrixElementPosition<T>>,
    tolerance: V) -> Result<SeparatedMatrix<T, V>, &'a str>
    where T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Copy + Debug + Eq + Hash + SubAssign + PartialOrd + AddAssign +
             From<u8> + 'static,
          V: Add<Output = V> + Mul<Output = V> + Sub<Output = V> + Div<Output = V> + Copy + Debug +
             PartialEq + AddAssign + MulAssign + SubAssign + Into<f64> + From<f32> + 'static
{
    let shape = matrix.get_shape();

    let all_elements_values =
        matrix.extract_all_elements_values();

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
            if positions.iter().position(|p| p.row() == i).is_none() &&
                positions.iter().position(|p| p.column() == j).is_none()
            {
                let value = extract_element_value(i, j, &all_elements_values);
                k_aa_elements.push(value);
            }
            j += T::from(1u8);
        }
        i += T::from(1u8);
    }

    // for i in 0..shape.0.into()
    // {
    //     for j in 0..shape.1.into()
    //     {
    //         if positions.iter().position(|p|
    //                 p.row == T::from(i)).is_none() &&
    //             positions.iter().position(|p|
    //                 p.column == T::from(j)).is_none()
    //         {
    //             let row = T::from(i);
    //             let column = T::from(j);
    //             let value = extract_element_value(row, column, &all_elements_values);
    //             k_aa_elements.push(value);
    //         }
    //     }
    // }

    let k_aa_matrix = ExtendedMatrix::create(k_aa_rows_number,
        k_aa_columns_number, k_aa_elements, tolerance);

    let k_ab_rows_number = shape.0 - converted_positions_length;

    let k_ab_columns_number = converted_positions_length;

    let mut k_ab_elements = Vec::new();

    let mut i = T::from(0u8);
    while i < shape.0
    {
        if positions.iter().position(|p| p.row() == i).is_none()
        {
            for j in 0..positions.len()
            {
                let row = i;
                let column = positions[j].column();
                if column > shape.1
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kab \
                        could not be composed!");
                }
                let value = extract_element_value(row, column, &all_elements_values);
                k_ab_elements.push(value);
            }
        }
        i += T::from(1u8);
    }

    // for i in 0..shape.0.into()
    // {
    //     if positions.iter().position(|p| p.row == T::from(i)).is_none()
    //     {
    //         for j in 0..positions.len()
    //         {
    //             let row = T::from(i);
    //             let column = positions[j].column;
    //             if column > shape.1
    //             {
    //                 return Err("Extended matrix: Matrix could not be separated! Matrix Kab \
    //                     could not be composed!");
    //             }
    //             let value = extract_element_value(row, column, &all_elements_values);
    //             k_ab_elements.push(value);
    //         }
    //     }
    // }

    let k_ab_matrix = ExtendedMatrix::create(k_ab_rows_number,
        k_ab_columns_number, k_ab_elements, tolerance);

    let k_ba_rows_number = converted_positions_length;

    let k_ba_columns_number = shape.1 - converted_positions_length;

    let mut k_ba_elements = Vec::new();


    for i in 0..positions.len()
    {
        let mut j = T::from(0u8);
        while j < shape.1
        {
            if positions.iter().position(|p| p.column() == j).is_none()
            {
                let row = positions[i].row();
                let column = j;
                if row > shape.0
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kba \
                        could not be composed!");
                }
                let value = extract_element_value(row, column, &all_elements_values);
                k_ba_elements.push(value);
            }
            j += T::from(1u8);
        }
    }

    //  for i in 0..positions.len()
    // {
    //     for j in 0..shape.1.into()
    //     {
    //         if positions.iter().position(|p|
    //             p.column == T::from(j)).is_none()
    //         {
    //             let row = positions[i].row;
    //             let column = T::from(j);
    //             if row > shape.0
    //             {
    //                 return Err("Extended matrix: Matrix could not be separated! Matrix Kba \
    //                     could not be composed!");
    //             }
    //             let value = extract_element_value(row, column, &all_elements_values);
    //             k_ba_elements.push(value);
    //         }
    //     }
    // }

    let k_ba_matrix = ExtendedMatrix::create(k_ba_rows_number,
        k_ba_columns_number, k_ba_elements, tolerance);

    let k_bb_rows_number = converted_positions_length;

    let k_bb_columns_number = converted_positions_length;

    let mut k_bb_elements = Vec::new();

    for i in 0..positions.len()
    {
        for j in 0..positions.len()
        {
            let row = positions[i].row();
            let column = positions[j].column();
            if row > shape.0 || column > shape.1
            {
                return Err("Extended matrix: Matrix could not be separated! Matrix Kbb could \
                    not be composed!");
            }
            let value = extract_element_value(row, column, &all_elements_values);
            k_bb_elements.push(value);
        }
    }
    let k_bb_matrix = ExtendedMatrix::create(k_bb_rows_number,
        k_bb_columns_number, k_bb_elements, tolerance);
    Ok(SeparatedMatrix {
        k_aa: k_aa_matrix, k_ab: k_ab_matrix, k_ba: k_ba_matrix, k_bb: k_bb_matrix
    })
}
