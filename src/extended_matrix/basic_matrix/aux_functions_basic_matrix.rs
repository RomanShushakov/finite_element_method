use crate::extended_matrix::{BasicMatrix};
use crate::extended_matrix::{SymmetricMatrix, NonSymmetricMatrix};

use std::fmt::Debug;


pub fn matrix_size_check<'a, T>(inputted_row: T, inputted_column: T, matrix_size: (T, T))
    -> Result<(), &'a str>
    where T: PartialOrd + Copy
{
    if inputted_row >= matrix_size.0 || inputted_column >= matrix_size.1
    {
        return Err("Basic matrix: Inputted indexes are out of matrix size!");
    }
    Ok(())
}


pub fn extract_value_by_index<T, V>(requested_index: T, indexes: &[T], values: &[V]) -> V
    where T: Copy + PartialEq,
          V: Copy + Default
{
    if let Some(position) = indexes
        .iter().position(|index| *index == requested_index)
    {
        values[position]
    }
    else
    {
        Default::default()
    }
}


pub fn return_symmetric_matrix_struct<T, V>(boxed_struct: Box<dyn BasicMatrix<T, V>>)
    -> SymmetricMatrix<T, V>
    where T: Copy + 'static,
          V: Copy + 'static
{
    let matrix: &SymmetricMatrix<T, V> = match boxed_struct
        .as_any()
        .downcast_ref::<SymmetricMatrix<T, V>>()
        {
            Some(matrix) => matrix,
            None => panic!("Basic matrix: Matrix is not symmetric!!!"),
        };
    let rows_and_columns_number = matrix.rows_and_columns_number;
    let elements_indexes = matrix.elements_indexes.clone();
    let elements_values = matrix.elements_values.clone();
    SymmetricMatrix { rows_and_columns_number, elements_indexes, elements_values }
}


pub fn return_non_symmetric_matrix_struct<T, V>(boxed_struct: Box<dyn BasicMatrix<T, V>>)
    -> NonSymmetricMatrix<T, V>
    where T: Copy + Debug + Default + 'static,
          V: Copy + Debug + Default + 'static
{
    let matrix: &NonSymmetricMatrix<T, V> = match boxed_struct
        .as_any()
        .downcast_ref::<NonSymmetricMatrix<T, V>>()
        {
            Some(matrix) => matrix,
            None => panic!("Basic matrix: Matrix is symmetric!!!"),
        };
    let rows_number = matrix.rows_number;
    let columns_number = matrix.columns_number;
    let elements_indexes = matrix.elements_indexes.clone();
    let elements_values = matrix.elements_values.clone();
    NonSymmetricMatrix { rows_number, columns_number, elements_indexes, elements_values }
}
