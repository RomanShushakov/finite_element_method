use crate::extended_matrix::BasicMatrixTrait;
use crate::extended_matrix::
    {
        MatrixElementPosition, Shape, ZerosRowColumn, NonSymmetricMatrix
    };
use crate::extended_matrix::BasicMatrixType;
use crate::extended_matrix::{matrix_size_check, extract_value_by_index};


use crate::ElementsNumbers;

use std::ops::{Sub, Add, Mul, MulAssign, Div, Rem, SubAssign};
use std::any::Any;
use std::fmt::Debug;
use std::collections::{HashMap};
use std::hash::Hash;


#[derive(Debug, Clone)]
pub struct SymmetricMatrix<T, V>
{
    pub rows_and_columns_number: T,
    pub elements_indexes: Vec<T>,
    pub elements_values: Vec<V>,
}


impl<T, V> BasicMatrixTrait<T, V> for SymmetricMatrix<T, V>
    where T: Copy + PartialOrd + Sub<Output = T> + Add<Output = T> + Mul<Output = T> +
             Div<Output = T> + Debug + Rem<Output = T> + Eq + Hash + Into<ElementsNumbers> +
             From<ElementsNumbers> + SubAssign + Default + 'static,
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


    fn into_symmetric(self) -> Box<dyn BasicMatrixTrait<T, V>>
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
        let mut can_continue = true;
        while can_continue
        {
            if let Some(row_column) = self.find_zeros_row_column()
            {
                let zeros_row_column = ZerosRowColumn { row: row_column, column: row_column };
                zeros_rows_columns.push(zeros_row_column);
                self.remove_row_column(row_column);
            }
            else
            {
                can_continue = false;
            }
        }
        zeros_rows_columns
    }


    fn remove_selected_row(&mut self, row: T) -> Box<dyn BasicMatrixTrait<T, V>>
    {
        let symmetric_matrix = self.clone();
        let mut non_symmetric_matrix = symmetric_matrix.non_symmetric();
        non_symmetric_matrix.remove_row(row);
        Box::new(non_symmetric_matrix)
    }


    fn remove_selected_column(&mut self, column: T) -> Box<dyn BasicMatrixTrait<T, V>>
    {
        let symmetric_matrix = self.clone();
        let mut non_symmetric_matrix = symmetric_matrix.non_symmetric();
        non_symmetric_matrix.remove_column(column);
        Box::new(non_symmetric_matrix)
    }
}


impl<T, V> SymmetricMatrix<T, V>
    where T: Copy + Debug + Into<ElementsNumbers> + From<ElementsNumbers> + PartialEq +
             Mul<Output = T> + Add<Output = T> + PartialOrd + SubAssign + Div<Output = T> +
             Sub<Output = T> + Rem<Output = T>,
          V: Copy + Default
{
    fn find_zeros_row_column(&self) -> Option<T>
    {
        let mut zeros_row_column = None;
        let find_index = |row, column| self.elements_indexes
            .iter()
            .position(|index|
                {
                    *index == T::from(row) * self.rows_and_columns_number + T::from(column)
                });
        for row_column in (0..self.rows_and_columns_number.into()).rev()
        {
            if (0..self.rows_and_columns_number.into()).rev().all(|row|
                find_index(row, row_column) == None) &&
                (0..self.rows_and_columns_number.into()).rev().all(|column|
                find_index(row_column, column) == None)
            {
                zeros_row_column = Some(T::from(row_column));
            }
            if zeros_row_column != None
            {
                break;
            }
        }
        zeros_row_column
    }


    fn remove_row_column(&mut self, row_column: T)
    {
        for index in self.elements_indexes.as_mut_slice()
        {
            if *index >= row_column * self.rows_and_columns_number
            {
                *index -= self.rows_and_columns_number;
            }
        }
        for index in self.elements_indexes.as_mut_slice()
        {
            if *index % self.rows_and_columns_number > row_column
            {
                *index -= *index / self.rows_and_columns_number + T::from(1);
            }
            else
            {
                *index -= *index / self.rows_and_columns_number;
            }
        }
        self.rows_and_columns_number -= T::from(1);
    }


    fn non_symmetric(&self) -> NonSymmetricMatrix<T, V>
    {
        let non_symmetric_rows_number = self.rows_and_columns_number;
        let non_symmetric_columns_number = self.rows_and_columns_number;
        let mut non_symmetric_indexes = Vec::new();
        let mut non_symmetric_values = Vec::new();
        for (index, value) in self.elements_indexes.iter().zip(self.elements_values.iter())
        {
            non_symmetric_indexes.push(*index);
            non_symmetric_values.push(*value);
            let current_row = *index / non_symmetric_columns_number;
            let current_column = *index % non_symmetric_columns_number;
            if current_row != current_column
            {
                let symmetric_index = current_column * non_symmetric_columns_number + current_row;
                non_symmetric_indexes.push(symmetric_index);
                non_symmetric_values.push(*value);
            }
        }
        NonSymmetricMatrix
        {
            rows_number: non_symmetric_rows_number, columns_number: non_symmetric_columns_number,
            elements_indexes: non_symmetric_indexes, elements_values: non_symmetric_values,
        }
    }
}
