use std::ops::{Div, Rem, Mul, Add, AddAssign, Sub, SubAssign, MulAssign};
use std::fmt::Debug;
use std::hash::Hash;

use extended_matrix::one::One;
use extended_matrix::basic_matrix::basic_matrix::MatrixElementPosition;
use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::extract_element_value;

use crate::fem::{StiffnessGroup, SeparatedMatrix};
use crate::fem::{StiffnessType};
use crate::fem::{STIFFNESS_TYPES_NUMBER, GLOBAL_DOF};
use crate::ElementsNumbers;

use crate::TOLERANCE;


pub fn global_dof<T>() -> T
    where T: One + Default
{
    let mut global_dof = T::default();
    (0..GLOBAL_DOF).for_each(|_| global_dof += T::one());
    global_dof
}


pub fn compose_stiffness_sub_groups<'a, T>(global_group_position: T,
    global_group_columns_number: T, global_number_1: T, global_number_2: T)
    -> Result<Vec<StiffnessGroup<T>>, &'a str>
    where T: Copy + Debug + Div<Output = T> + Rem<Output = T> + Mul<Output = T> + Add<Output = T> +
             PartialOrd + AddAssign + Default
{
    let mut stiffness_sub_groups = Vec::new();
    let row = global_group_position / global_group_columns_number;
    let column = global_group_position % global_group_columns_number;
    let mut k = T::default();
    while k < T::from(STIFFNESS_TYPES_NUMBER as usize)
    {
        let start_row = row * T::from(GLOBAL_DOF as usize);
        let row_shift_init = k / T::from(2) * (T::from(GLOBAL_DOF as usize) / T::from(2));
        let row_shift_final = k / T::from(2) * (T::from(GLOBAL_DOF as usize) / T::from(2)) +
            (T::from(GLOBAL_DOF as usize) / T::from(2));
        let start_column = column * T::from(GLOBAL_DOF as usize);
        let column_shift_init = k % T::from(2) * (T::from(GLOBAL_DOF as usize) / T::from(2));
        let column_shift_final = k % T::from(2) * (T::from(GLOBAL_DOF as usize) / T::from(2)) +
            (T::from(GLOBAL_DOF as usize) / T::from(2));
        let mut element_positions = Vec::new();
        let mut current_row = start_row + row_shift_init;
        while current_row < start_row + row_shift_final
        {
            let mut current_column  = start_column + column_shift_init;
            while current_column < start_column + column_shift_final
            {
                element_positions.push(MatrixElementPosition { row: current_row,
                    column: current_column });
                current_column += T::from(1);
            }
            current_row += T::from(1);
        }
        let stiffness_type = StiffnessType::iterator()
            .nth(k.into() as usize)
            .ok_or("FEModel: Stiffness type could not be defined")?;
        let stiffness_sub_group = StiffnessGroup { stiffness_type: *stiffness_type,
            number_1: global_number_1,
            number_2: global_number_2,
            positions: element_positions,
        };
        stiffness_sub_groups.push(stiffness_sub_group);
        k += T::from(1);
    }
    Ok(stiffness_sub_groups)
}


pub fn separate<'a, T, V>(matrix: ExtendedMatrix<T, V>, positions: Vec<MatrixElementPosition<T>>)
    -> Result<SeparatedMatrix<T, V>, &'a str>
    where T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Copy + Default + Debug + Eq + Hash + SubAssign + PartialOrd +
             'static,
          V: Add<Output = V> + Mul<Output = V> + Sub<Output = V> + Div<Output = V> + Copy + Debug +
             PartialEq + Default + AddAssign + MulAssign + SubAssign + Into<f64> + 'static
{
    let shape = matrix.get_shape();
    let all_elements_values =
        matrix.extract_all_elements_values();
    let k_aa_rows_number = shape.0 - T::from(positions.len());
    let k_aa_columns_number = shape.1 - T::from(positions.len());
    let mut k_aa_elements = Vec::new();
    for i in 0..shape.0.into()
    {
        for j in 0..shape.1.into()
        {
            if positions.iter().position(|p|
                    p.row == T::from(i)).is_none() &&
                positions.iter().position(|p|
                    p.column == T::from(j)).is_none()
            {
                let row = T::from(i);
                let column = T::from(j);
                let value = extract_element_value(row, column, &all_elements_values);
                k_aa_elements.push(value);
            }
        }
    }
    let k_aa_matrix =
        ExtendedMatrix::create(k_aa_rows_number, k_aa_columns_number, k_aa_elements,
            TOLERANCE);
    let k_ab_rows_number = shape.0 - T::from(positions.len());
    let k_ab_columns_number = T::from(positions.len());
    let mut k_ab_elements = Vec::new();
    for i in 0..shape.0.into()
    {
        if positions.iter().position(|p| p.row == T::from(i)).is_none()
        {
            for j in 0..positions.len()
            {
                let row = T::from(i);
                let column = positions[j].column;
                if column > shape.1
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kab \
                        could not be composed!");
                }
                let value = extract_element_value(row, column, &all_elements_values);
                k_ab_elements.push(value);
            }
        }
    }
    let k_ab_matrix =
        ExtendedMatrix::create(k_ab_rows_number, k_ab_columns_number, k_ab_elements,
            TOLERANCE);
    let k_ba_rows_number = T::from(positions.len());
    let k_ba_columns_number = shape.1 - T::from(positions.len());
    let mut k_ba_elements = Vec::new();
    for i in 0..positions.len()
    {
        for j in 0..shape.1.into()
        {
            if positions.iter().position(|p|
                p.column == T::from(j)).is_none()
            {
                let row = positions[i].row;
                let column = T::from(j);
                if row > shape.0
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kba \
                        could not be composed!");
                }
                let value = extract_element_value(row, column, &all_elements_values);
                k_ba_elements.push(value);
            }
        }
    }
    let k_ba_matrix =
        ExtendedMatrix::create(k_ba_rows_number, k_ba_columns_number, k_ba_elements,
            TOLERANCE);
    let k_bb_rows_number = T::from(positions.len());
    let k_bb_columns_number = T::from(positions.len());
    let mut k_bb_elements = Vec::new();
    for i in 0..positions.len()
    {
        for j in 0..positions.len()
        {
            let row = positions[i].row;
            let column = positions[j].column;
            if row > shape.0 || column > shape.1
            {
                return Err("Extended matrix: Matrix could not be separated! Matrix Kbb could \
                    not be composed!");
            }
            let value = extract_element_value(row, column, &all_elements_values);
            k_bb_elements.push(value);
        }
    }
    let k_bb_matrix =
        ExtendedMatrix::create(k_bb_rows_number, k_bb_columns_number, k_bb_elements,
            TOLERANCE);
    Ok(SeparatedMatrix {
        k_aa: k_aa_matrix, k_ab: k_ab_matrix, k_ba: k_ba_matrix, k_bb: k_bb_matrix
    })
}
