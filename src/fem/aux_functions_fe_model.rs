use crate::fem::{StiffnessGroup};
use crate::fem::{StiffnessType};
use crate::fem::{STIFFNESS_TYPES_NUMBER, GLOBAL_DOF};
use crate::ElementsNumbers;

use crate::extended_matrix::MatrixElementPosition;

use std::ops::{Div, Rem, Mul, Add, AddAssign};
use std::fmt::Debug;


pub fn compose_stiffness_sub_groups<'a, T>(global_group_position: T,
    global_group_columns_number: T, global_number_1: T, global_number_2: T)
    -> Result<Vec<StiffnessGroup<T>>, &'a str>
    where T: Copy + Debug + Div<Output = T> + Rem<Output = T> + Mul<Output = T> +
             From<ElementsNumbers> + Into<ElementsNumbers> + Add<Output = T> + PartialOrd +
             AddAssign + Default
{
    let mut stiffness_sub_groups = Vec::new();
    let row = global_group_position / global_group_columns_number;
    let column = global_group_position % global_group_columns_number;
    let mut k = T::default();
    while k < STIFFNESS_TYPES_NUMBER.into()
    {
        let start_row = row * T::from(GLOBAL_DOF);
        let row_shift_init = k / T::from(2) * (T::from(GLOBAL_DOF) / T::from(2));
        let row_shift_final = k / T::from(2) * (T::from(GLOBAL_DOF) / T::from(2)) +
            (T::from(GLOBAL_DOF) / T::from(2));
        let start_column = column * T::from(GLOBAL_DOF);
        let column_shift_init = k % T::from(2) * (T::from(GLOBAL_DOF) / T::from(2));
        let column_shift_final = k % T::from(2) * (T::from(GLOBAL_DOF) / T::from(2)) +
            (T::from(GLOBAL_DOF) / T::from(2));
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