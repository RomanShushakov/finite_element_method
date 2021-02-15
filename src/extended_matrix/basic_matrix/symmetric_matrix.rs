use crate::extended_matrix::basic_matrix::{BasicMatrix};
use crate::extended_matrix::basic_matrix::{Shape, MatrixElementPosition, ZerosRowColumn};
use crate::extended_matrix::basic_matrix::{BasicMatrixType};
use crate::extended_matrix::basic_matrix::{matrix_size_check, extract_value_by_index};
use crate::ElementsNumbers;

use std::ops::{Sub, Add, Mul, MulAssign, Div, Rem};
use std::any::Any;
use std::fmt::Debug;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;


#[derive(Debug, Clone)]
pub struct SymmetricMatrix<T, V>
{
    pub rows_and_columns_number: T,
    pub elements_indexes: Vec<T>,
    pub elements_values: Vec<V>,
}


impl<T, V> BasicMatrix<T, V> for SymmetricMatrix<T, V>
    where T: Copy + PartialOrd + Sub<Output = T> + Add<Output = T> + Mul<Output = T> +
             Div<Output = T> + Debug + Rem<Output = T> + Eq + Hash + Into<ElementsNumbers> +
             From<ElementsNumbers> + 'static,
          V: Copy + Default + Debug + PartialEq + MulAssign + 'static,
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
            (self.rows_and_columns_number, self.rows_and_columns_number))?;
        let (row, column) = if row <= column { (row, column) } else { (column, row) };
        let requested_index = row * self.rows_and_columns_number + column;
        let value = extract_value_by_index(
            requested_index, self.elements_indexes.as_slice(), self.elements_values.as_slice());
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
    //         (self.rows_and_columns_number, self.rows_and_columns_number))?;
    //     let (row, column) = if row <= column { (row, column) } else { (column, row) };
    //     let requested_index = row * self.rows_and_columns_number + column;
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
    //         (self.rows_and_columns_number, self.rows_and_columns_number))?;
    //     let (row, column) = if row <= column { (row, column) } else { (column, row) };
    //     let requested_index = row * self.rows_and_columns_number + column;
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
            let row = *index / self.rows_and_columns_number;
            let column = *index % self.rows_and_columns_number;
            let position = MatrixElementPosition { row, column };
            all_elements_values.insert(position, *value);
            if row != column
            {
                let symmetric_position = MatrixElementPosition { row: column, column: row };
                all_elements_values.insert(symmetric_position, *value);
            }
        }
        all_elements_values
    }


    fn get_shape(&self) -> Shape<T>
    {
        Shape(self.rows_and_columns_number, self.rows_and_columns_number)
    }


    fn transpose(&mut self) { }


    fn multiply_by_number(&mut self, number: V)
    {
        for i in 0..self.elements_values.len()
        {
            self.elements_values[i] *= number;
        }
    }


    fn into_symmetric(self) -> Box<dyn BasicMatrix<T, V>>
    {
        Box::new(self)
    }



    fn define_type(&self) -> BasicMatrixType
    {
        BasicMatrixType::Symmetric
    }


    fn as_any(&self) -> &dyn Any
    {
        self
    }


    fn remove_zeros_rows_columns(&mut self) -> Vec<ZerosRowColumn<T>>
    {
        let mut zeros_rows_columns = Vec::new();
        let mut non_zeros_rows = HashSet::new();
        let mut non_zeros_columns = HashSet::new();
        for index in &self.elements_indexes
        {
            let non_zeros_row = *index / self.rows_and_columns_number;
            let non_zeros_column = *index % self.rows_and_columns_number;
            non_zeros_rows.insert(non_zeros_row);
            non_zeros_columns.insert(non_zeros_column);
        }
        for i in 0..self.rows_and_columns_number.into()
        {
            for j in 0..self.rows_and_columns_number.into()
            {
                if non_zeros_rows.get(&T::from(i)) == None &&
                    non_zeros_columns.get(&T::from(j)) == None
                {
                    let zeros_row_column = ZerosRowColumn { row: T::from(i), column: T::from(j) };
                    zeros_rows_columns.push(zeros_row_column);
                    if i != j
                    {
                        let symmetric_zeros_row_column = ZerosRowColumn { row: T::from(j), column: T::from(i) };
                        zeros_rows_columns.push(symmetric_zeros_row_column);
                    }
                }
            }
        }
        zeros_rows_columns
    }
}
