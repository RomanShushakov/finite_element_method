use crate::extended_matrix::basic_matrix::{BasicMatrix};
use crate::extended_matrix::basic_matrix::{SymmetricMatrix, Shape, MatrixElementPosition};
use crate::extended_matrix::basic_matrix::{BasicMatrixType};
use crate::extended_matrix::basic_matrix::{matrix_size_check, extract_value_by_index};

use std::ops::{Mul, Add, Sub, Div, Rem, MulAssign};
use std::fmt::Debug;
use std::any::Any;
use std::collections::HashMap;
use std::hash::Hash;


#[derive(Debug, Clone)]
pub struct NonSymmetricMatrix<T, V>
    where T: Copy + Debug,
          V: Copy + Default
{
    pub rows_number: T,
    pub columns_number: T,
    pub elements_indexes: Vec<T>,
    pub elements_values: Vec<V>,
}


impl<T, V> BasicMatrix<T, V> for NonSymmetricMatrix<T, V>
    where T: Copy + PartialEq + Debug + PartialOrd + Mul<Output = T> +
             Add<Output = T> + Default + Sub<Output = T> + Div<Output = T> + Rem<Output = T> +
             Eq + Hash + 'static,
          V: Copy + Default + PartialEq + Debug + MulAssign + 'static,
{
   // fn create_element_value(&mut self, requested_index: T, new_value: V)
    // {
    //     self.elements_indexes.push(requested_index);
    //     self.elements_values.push(new_value);
    // }


    fn read_element_value(&self, row: T, column: T) -> Result<V, &str>
    {
        matrix_size_check(
            row, column,
            (self.rows_number, self.columns_number))?;
        let requested_index = row * self.columns_number + column;
        let value = extract_value_by_index(
            requested_index,
            self.elements_indexes.as_slice(),
            self.elements_values.as_slice());
        Ok(value)
    }


    // fn update_element_value(&mut self, row: T, column: T, new_value: V) -> Result<(), &str>
    // {
    //     if new_value == Default::default()
    //     {
    //         self.delete_element_value(row, column)?;
    //         return Ok(());
    //     }
    //     matrix_size_check(
    //         row, column,
    //         (self.rows_number, self.columns_number))?;
    //     let requested_index = row * self.columns_number + column;
    //     if let Some(position) = self.elements_indexes
    //         .iter().position(|index| *index == requested_index)
    //     {
    //         self.elements_values[position] = new_value;
    //     }
    //     else
    //     {
    //         self.create_element_value(requested_index, new_value);
    //     }
    //     Ok(())
    // }


    // fn delete_element_value(&mut self, row: T, column: T) -> Result<(), &str>
    // {
    //     matrix_size_check(
    //         row, column,
    //         (self.rows_number, self.columns_number))?;
    //     let requested_index = row * self.columns_number + column;
    //     if let Some(position) = self.elements_indexes
    //         .iter().position(|index| *index == requested_index)
    //     {
    //         self.elements_indexes.remove(position);
    //         self.elements_values.remove(position);
    //     }
    //     Ok(())
    // }


    fn extract_all_elements_values(&self) -> HashMap<MatrixElementPosition<T>, V>
    {
        let mut all_elements_values = HashMap::new();
        for (index, value) in self.elements_indexes.iter()
            .zip(self.elements_values.iter())
        {
            let row = *index / self.columns_number;
            let column = *index % self.columns_number;
            let position = MatrixElementPosition { row, column };
            all_elements_values.insert(position, *value);
        }
        all_elements_values
    }


    fn get_shape(&self) -> Shape<T>
    {
        Shape(self.rows_number, self.columns_number)
    }


    fn transpose(&mut self)
    {
        let transposed_rows_number = self.columns_number;
        let transposed_columns_number = self.rows_number;
        for i in 0..self.elements_indexes.len()
        {
            let current_index_row = self.elements_indexes[i] / self.columns_number;
            let current_index_column = self.elements_indexes[i] % self.columns_number;
            let transported_index_row = current_index_column;
            let transported_index_column = current_index_row;
            let transported_index = transported_index_row * transposed_columns_number +
                transported_index_column;
            self.elements_indexes[i] = transported_index;
        }
        self.rows_number = transposed_rows_number;
        self.columns_number = transposed_columns_number;
    }


    fn multiply_by_number(&mut self, number: V)
    {
        for i in 0..self.elements_values.len()
        {
            self.elements_values[i] *= number;
        }
    }


    fn into_symmetric(self) -> Box<dyn BasicMatrix<T, V>>
    {
        if self.rows_number != self.columns_number
        {
            return Box::new(self);
        }
        let mut elements_indexes = Vec::new();
        let mut elements_values = Vec::new();
        let mut indexes = self.elements_indexes.clone();
        let mut values = self.elements_values.clone();
        while !indexes.is_empty()
        {
            let current_index = indexes.remove(0);
            let current_value = values.remove(0);
            let current_row = current_index / self.rows_number;
            let current_column = current_index % self.rows_number;
            if current_row == current_column
            {
                elements_indexes.push(current_index);
                elements_values.push(current_value);
                continue;
            }
            if let Some(position) = indexes
                .iter()
                .position(|index| *index == current_column * self.rows_number + current_row)
            {
                let current_symmetric_value = values.remove(position);
                if current_value != current_symmetric_value
                {
                    return Box::new(self);
                }
                indexes.remove(position);
                let (row, column) =
                    if current_row < current_column { (current_row, current_column) }
                    else { (current_column, current_row) };
                elements_indexes.push(row * self.rows_number + column);
                elements_values.push(current_value);
            }
            else
            {
                return Box::new(self);
            }
        }
        Box::new(SymmetricMatrix
            {
                rows_and_columns_number: self.rows_number, elements_indexes, elements_values
            })
    }


    fn define_type(&self) -> BasicMatrixType
    {
        BasicMatrixType::NonSymmetric
    }


    fn as_any(&self) -> &dyn Any
    {
        self
    }
}
