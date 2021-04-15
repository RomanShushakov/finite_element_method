use crate::extended_matrix::BasicMatrixTrait;
use crate::extended_matrix::{MatrixElementPosition, ZerosRowColumn, Shape, SymmetricMatrix};
use crate::extended_matrix::BasicMatrixType;
use crate::extended_matrix::
    {
        matrix_size_check, extract_value_by_index
    };

use crate::ElementsNumbers;

use std::ops::{Mul, Add, Sub, Div, Rem, MulAssign, SubAssign};
use std::fmt::Debug;
use std::any::Any;
use std::collections::{HashMap};
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


impl<T, V> BasicMatrixTrait<T, V> for NonSymmetricMatrix<T, V>
    where T: Copy + PartialEq + Debug + PartialOrd + Mul<Output = T> + Add<Output = T> +
             Default + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash +
             Into<ElementsNumbers> + From<ElementsNumbers> + SubAssign + 'static,
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


    fn into_symmetric(self) -> Box<dyn BasicMatrixTrait<T, V>>
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


    fn remove_zeros_rows_columns(&mut self) -> Vec<ZerosRowColumn<T>>
    {
        let mut zeros_rows_columns = Vec::new();
        let mut can_continue = true;
        while can_continue
        {
            if let Some(row) = self.find_zeros_row()
            {
                if let Some(column) = self.find_zeros_column()
                {
                    let zeros_row_column = ZerosRowColumn { row, column };
                    zeros_rows_columns.push(zeros_row_column);
                    self.remove_row(row);
                    self.remove_column(column);
                }
                else
                {
                    can_continue = false;
                }
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
        self.remove_row(row);
        Box::new(self.clone())
    }


    fn remove_selected_column(&mut self, column: T) -> Box<dyn BasicMatrixTrait<T, V>>
    {
        self.remove_column(column);
        Box::new(self.clone())
    }
}


impl<T, V> NonSymmetricMatrix<T, V>
    where T: Copy + Debug + Into<ElementsNumbers> + From<ElementsNumbers> + PartialEq +
             Mul<Output = T> + Add<Output = T> + PartialOrd + SubAssign + Div<Output = T> +
             Sub<Output = T> + Rem<Output = T>,
          V: Copy + Default
{
    fn find_zeros_row(&self) -> Option<T>
    {
        let mut zeros_row = None;
        let find_index = |row, column| self.elements_indexes
            .iter()
            .position(|index|
                {
                    *index == T::from(row) * self.columns_number + T::from(column)
                });
        for row in (0..self.rows_number.into()).rev()
        {
            if (0..self.columns_number.into()).rev().all(|column|
                find_index(row, column) == None)
            {
                zeros_row = Some(T::from(row));
            }
            if zeros_row != None
            {
                break;
            }
        }
        zeros_row
    }


    fn find_zeros_column(&self) -> Option<T>
    {
        let mut zeros_column = None;
        let find_index = |row, column| self.elements_indexes
            .iter()
            .position(|index| *index == T::from(row) * self.columns_number + T::from(column));
        for column in (0..self.columns_number.into()).rev()
        {
            if (0..self.rows_number.into()).rev().all(|row| find_index(row, column) == None)
            {
                zeros_column = Some(T::from(column));
            }
            if zeros_column != None
            {
                break;
            }
        }
        zeros_column
    }


    pub fn remove_row(&mut self, row: T)
    {
        for index in self.elements_indexes.as_mut_slice()
        {
            if *index > row * self.columns_number
            {
                *index -= self.columns_number;
            }
        }
        self.rows_number -= T::from(1);
    }


    pub fn remove_column(&mut self, column: T)
    {
        for index in self.elements_indexes.as_mut_slice()
        {
            if *index % self.columns_number > column
            {
                *index -= *index / self.columns_number + T::from(1);
            }
            else
            {
                *index -= *index / self.columns_number;
            }
        }
        self.columns_number -= T::from(1);
    }
}
