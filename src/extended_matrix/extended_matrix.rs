use crate::extended_matrix::BasicMatrixTrait;
use crate::extended_matrix::
    {
        NonSymmetricMatrix, MatrixElementPosition, ZerosRowColumn
    };
use crate::extended_matrix::
    {
        matrices_dimensions_conformity_check, extract_element_value, remove_zero_values
    };

use crate::fem::SeparatedMatrix;

use crate::{ElementsNumbers, ElementsValues, TOLERANCE};

use std::fmt::Debug;
use std::convert::{From, Into};
use std::ops::{Mul, Add, Sub, Div, Rem, MulAssign, AddAssign, SubAssign};
use std::hash::Hash;


#[derive(Copy, Clone)]
pub enum Operation
{
    Addition,
    Multiplication,
    Subtraction,
}


#[derive(Clone)]
pub struct ExtendedMatrix<T, V>
{
    pub basic_matrix: Box<dyn BasicMatrixTrait<T, V>>
}


impl<T, V> ExtendedMatrix<T, V>
    where T: Copy + From<ElementsNumbers> + Into<ElementsNumbers> + Debug + Mul<Output = T> +
             PartialOrd + Add + Sub + Add<Output = T> + Sub<Output = T> + Default + Div + Rem +
             Div<Output = T> + Rem<Output = T> + Eq + Hash + SubAssign + 'static,
          V: Copy + Debug + PartialEq + Default + AddAssign + MulAssign + Mul<Output = V> +
             Div<Output = V> + SubAssign + Sub<Output = V> + Into<ElementsValues> +
             From<ElementsValues> + Add<Output = V> + 'static,
{
    pub fn show_matrix(&self)
    {
        let shape = self.basic_matrix.get_shape();
        for i in 0..shape.0.into()
        {
            for j in 0..shape.1.into()
            {
                print!("{:?}, ", self.basic_matrix
                    .read_element_value(T::from(i), T::from(j)).unwrap());
            }
            println!();
        }
    }


    pub fn create(rows_number: T, columns_number: T, all_elements: Vec<V>) -> Self
    {
        let mut elements_indexes = Vec::new();
        let mut elements_values = Vec::new();
        let all_elements_length = rows_number * columns_number;
        for (index, value) in (0..all_elements_length.into()).zip(all_elements)
        {
            if value.into().abs() > TOLERANCE
            {
                elements_indexes.push(T::from(index));
                elements_values.push(value);
            }
        }
        let basic_matrix = Box::new(NonSymmetricMatrix
            {
                columns_number, rows_number, elements_indexes, elements_values
            });
        let basic_matrix = basic_matrix.into_symmetric();
        ExtendedMatrix { basic_matrix }
    }


    pub fn transpose(&mut self)
    {
        self.basic_matrix.transpose();
    }


    pub fn add_subtract_matrix<'a>(&'a self, other: &'a Self, operation: Operation)
        -> Result<Self, &'a str>
    {
        let (_, shape) =
            matrices_dimensions_conformity_check(&self, &other, operation)?;
        let lhs_all_elements_values =
            self.basic_matrix.extract_all_elements_values();
        let rhs_all_elements_values =
            other.basic_matrix.extract_all_elements_values();
        let mut elements_indexes = Vec::new();
        let mut elements_values = Vec::new();
        for index in 0..(shape.0 * shape.1).into()
        {
            let current_lhs_element_value = extract_element_value(
                    T::from(index) / shape.1, T::from(index) % shape.1,
                    &lhs_all_elements_values
                );
            let current_rhs_element_value = extract_element_value(
                    T::from(index) / shape.1, T::from(index) % shape.1,
                    &rhs_all_elements_values
                );
            let value =
                {
                    match operation
                    {
                        Operation::Addition =>
                            current_lhs_element_value + current_rhs_element_value,
                        Operation::Subtraction =>
                            current_lhs_element_value - current_rhs_element_value,
                        Operation::Multiplication =>
                            return Err("Extended matrix: Multiplication operation could not be \
                                applied for add_subtract function!"),
                    }
                };
            if value.into().abs() > TOLERANCE
            {
                elements_indexes.push(T::from(index));
                elements_values.push(value);
            }
        }
        let basic_matrix = Box::new(NonSymmetricMatrix
            {
                rows_number: shape.0, columns_number: shape.1, elements_indexes, elements_values
            });
        let basic_matrix = basic_matrix.into_symmetric();
        Ok(ExtendedMatrix { basic_matrix })
    }


    pub fn add_matrix<'a>(&'a self, other: &'a Self) -> Result<Self, &'a str>
    {
        self.add_subtract_matrix(other, Operation::Addition)
    }


    pub fn subtract_matrix<'a>(&'a self, other: &'a Self) -> Result<Self, &'a str>
    {
        self.add_subtract_matrix(other, Operation::Subtraction)
    }


    pub fn add_sub_matrix(&mut self, other: &Self, self_positions: &[MatrixElementPosition<T>],
        other_positions: &[MatrixElementPosition<T>])
    {
        let lhs_shape = self.basic_matrix.get_shape();
        let lhs_all_elements_values =
            self.basic_matrix.extract_all_elements_values();
        let rhs_all_elements_values =
            other.basic_matrix.extract_all_elements_values();
        let mut elements_indexes = Vec::new();
        let mut elements_values = Vec::new();
        for (lhs_position, rhs_position) in
            self_positions.iter().zip(other_positions)
        {
            let current_lhs_element_value = extract_element_value(
                    lhs_position.row, lhs_position.column,
                    &lhs_all_elements_values
                );
            let current_rhs_element_value = extract_element_value(
                    rhs_position.row, rhs_position.column,
                    &rhs_all_elements_values
                );
            let value = current_lhs_element_value + current_rhs_element_value;
            if value.into().abs() > TOLERANCE
            {
                elements_indexes.push(lhs_position.row * lhs_shape.1 + lhs_position.column);
                elements_values.push(value);
            }
        }
        for index in 0..(lhs_shape.0 * lhs_shape.1).into()
        {
            if let None = self_positions.iter().position(|position|
                *position == MatrixElementPosition { row: T::from(index) / lhs_shape.1,
                    column: T::from(index) % lhs_shape.1 })
            {
                let value = extract_element_value(
                    T::from(index) / lhs_shape.1, T::from(index) % lhs_shape.1,
                    &lhs_all_elements_values
                );
                if value.into().abs() > TOLERANCE
                {
                    elements_indexes.push(T::from(index));
                    elements_values.push(value);
                }
            }
        }
        let basic_matrix = Box::new(NonSymmetricMatrix
            {
                rows_number: lhs_shape.0, columns_number: lhs_shape.1, elements_indexes,
                elements_values
            });
        let basic_matrix = basic_matrix.into_symmetric();
        self.basic_matrix = basic_matrix;
    }


    pub fn multiply_by_number(&mut self, number: V)
    {
        self.basic_matrix.multiply_by_number(number);
    }


    pub fn multiply_by_matrix<'a>(&'a self, other: &'a Self) -> Result<Self, &'a str>
    {
        let (basic_dimension, shape) =
            matrices_dimensions_conformity_check(&self, &other,
        Operation::Multiplication)?;
        let lhs_all_elements_values =
            self.basic_matrix.extract_all_elements_values();
        let rhs_all_elements_values =
            other.basic_matrix.extract_all_elements_values();
        let mut elements_indexes = Vec::new();
        let mut elements_values = Vec::new();
        for index in 0..(shape.0 * shape.1).into()
        {
            let mut value = V::default();
            for k in 0..basic_dimension.into()
            {
                let current_lhs_element_value = extract_element_value(
                        T::from(index) / shape.1, T::from(k),
                        &lhs_all_elements_values
                    );
                let current_rhs_element_value = extract_element_value(
                        T::from(k), T::from(index) % shape.1,
                        &rhs_all_elements_values
                    );
                value += current_lhs_element_value * current_rhs_element_value;
            }
            if value.into().abs() > TOLERANCE
            {
                elements_indexes.push(T::from(index));
                elements_values.push(value);
            }
        }
        let basic_matrix = Box::new(NonSymmetricMatrix
            {
                rows_number: shape.0, columns_number: shape.1, elements_indexes, elements_values
            });
        let basic_matrix = basic_matrix.into_symmetric();
        Ok(ExtendedMatrix { basic_matrix })
    }


    pub fn naive_gauss_elimination<'a>(&'a self, other: &'a Self) -> Result<Self, &'a str>
    {
        let (basic_dimension, shape) =
            matrices_dimensions_conformity_check(&self, &other,
             Operation::Multiplication)?;
        let mut lhs_all_elements_values =
            self.basic_matrix.extract_all_elements_values();
        let mut rhs_all_elements_values =
            other.basic_matrix.extract_all_elements_values();
        let mut elements_values = (0..shape.0.into()).map(|_| V::default())
            .collect::<Vec<V>>();
        for k in 0..basic_dimension.into() - 1
        {
            for i in (k + 1)..basic_dimension.into()
            {
                let current_lhs_element_value = extract_element_value(
                        T::from(i), T::from(k), &lhs_all_elements_values
                    );
                let current_diag_lhs_element_value = extract_element_value(
                        T::from(k), T::from(k), &lhs_all_elements_values
                    );
                let current_coefficient = current_lhs_element_value / current_diag_lhs_element_value;
                for j in (k + 1)..basic_dimension.into()
                {
                    let current_lhs_element_value = extract_element_value(
                            T::from(k), T::from(j), &lhs_all_elements_values
                        );
                    *lhs_all_elements_values
                        .entry(MatrixElementPosition { row: T::from(i), column: T::from(j) })
                        .or_insert(Default::default()) -=
                        current_coefficient * current_lhs_element_value;
                }
                let current_rhs_element_value = extract_element_value(
                        T::from(k), T::default(), &rhs_all_elements_values
                    );
                *rhs_all_elements_values
                    .entry(MatrixElementPosition { row: T::from(i), column: T::default() })
                    .or_insert(Default::default()) -=
                    current_coefficient * current_rhs_element_value;
            }
        }
        let rhs_element_value = extract_element_value(
                T::from(basic_dimension.into() - 1), T::default(),
                &rhs_all_elements_values
            );
        let lhs_element_value = extract_element_value(
                T::from(basic_dimension.into() - 1), T::from(basic_dimension.into() - 1),
                &lhs_all_elements_values
            );
        elements_values[(basic_dimension.into() - 1) as usize] =
            rhs_element_value / lhs_element_value;
        for i in (0..basic_dimension.into() - 1).into_iter().rev()
        {
            let rhs_element_value = extract_element_value(
                    T::from(i), T::default(), &rhs_all_elements_values
                );
            let mut sum = rhs_element_value;
            for j in (i + 1)..basic_dimension.into()
            {
                let lhs_element_value = extract_element_value(
                        T::from(i), T::from(j), &lhs_all_elements_values
                    );
                sum -= lhs_element_value * elements_values[j as usize];
            }
            let lhs_element_value = extract_element_value(
                    T::from(i), T::from(i), &lhs_all_elements_values
                );
            elements_values[i as usize] = sum / lhs_element_value;
        }
        Ok(ExtendedMatrix::create(shape.0, shape.1, elements_values))
    }


    pub fn lu_decomposition(&self) -> Result<(Self, Self), &str>
    {
        let shape = self.basic_matrix.get_shape();
        if (shape.0 != shape.1) || shape.0 < T::from(2)
        {
            return Err("Extended matrix: Matrix could not be decomposed!");
        }
        let mut l_elements_indexes = Vec::new();
        let mut l_elements_values= Vec::new();
        for i in 0..shape.0.into()
        {
            l_elements_indexes.push(T::from(i) * shape.1 + T::from(i));
            l_elements_values.push(V::from(1.0));
        }
        let mut all_elements_values =
            self.basic_matrix.extract_all_elements_values();
        let mut u_elements_indexes = Vec::new();
        let mut u_elements_values= Vec::new();
        for k in 0..shape.1.into()
        {
             let current_element_value = extract_element_value(
                T::default(), T::from(k), &all_elements_values);
            u_elements_indexes.push(T::from(k));
            u_elements_values.push(current_element_value);
        }
        let mut row_number = 0;
        while row_number < shape.0.into() - 1
        {
            for i in (row_number + 1)..shape.0.into()
            {
                let current_coefficient = extract_element_value(
                    T::from(i), T::from(row_number), &all_elements_values) /
                    extract_element_value(T::from(row_number), T::from(row_number),
                    &all_elements_values);
                l_elements_indexes.push(T::from(i) * shape.1 + T::from(row_number));
                l_elements_values.push(current_coefficient);
                for j in 0..shape.1.into()
                {
                    let current_element_value = extract_element_value(
                        T::from(i), T::from(j), &all_elements_values) -
                        extract_element_value(T::from(row_number), T::from(j),
                    &all_elements_values) * current_coefficient;
                    if let Some(position) = u_elements_indexes
                        .iter()
                        .position(|index| *index ==  T::from(i) * shape.1 + T::from(j))
                    {
                        u_elements_values[position] = current_element_value;
                    }
                    else
                    {
                        u_elements_indexes.push(T::from(i) * shape.1 + T::from(j));
                        u_elements_values.push(current_element_value);
                    }
                    *all_elements_values
                        .entry(MatrixElementPosition {
                            row: T::from(i), column: T::from(j) })
                        .or_insert(Default::default()) = current_element_value;
                }
            }
            row_number += 1;
        }
        remove_zero_values(&mut l_elements_indexes, &mut l_elements_values);
        let l_basic_matrix = Box::new(NonSymmetricMatrix
            {
                rows_number: shape.0, columns_number: shape.1,
                elements_indexes: l_elements_indexes, elements_values: l_elements_values
            });
        let l_matrix = ExtendedMatrix { basic_matrix: l_basic_matrix };
        remove_zero_values(&mut u_elements_indexes, &mut u_elements_values);
        let u_basic_matrix = Box::new(NonSymmetricMatrix
            {
                rows_number: shape.0, columns_number: shape.1,
                elements_indexes: u_elements_indexes, elements_values: u_elements_values
            });
        let u_matrix = ExtendedMatrix { basic_matrix: u_basic_matrix };
        Ok((l_matrix, u_matrix))
    }


    pub fn determinant(&self) -> Result<V, &str>
    {
        let (_, u_matrix) = self.lu_decomposition()?;
        let u_matrix_elements_values = u_matrix.basic_matrix
            .extract_all_elements_values();
        let shape = u_matrix.basic_matrix.get_shape();
        let mut determinant = V::from(1.0);
        for i in 0..shape.0.into()
        {
            let current_diag_element_value = extract_element_value(
                T::from(i), T::from(i), &u_matrix_elements_values);
            determinant *= current_diag_element_value;
        }
        Ok(determinant)
    }


    pub fn inverse(&self) -> Result<Self, &str>
    {
        let (l_matrix, u_matrix) = self.lu_decomposition()?;
        let shape = self.basic_matrix.get_shape();
        let mut inverse_matrix_indexes = Vec::new();
        let mut inverse_matrix_values = Vec::new();
        for k in 0..shape.1.into()
        {
            let unit_column_indexes = vec![T::from(k)];
            let unit_column_values = vec![V::from(1.0)];
            let basic_unit_column = Box::new(NonSymmetricMatrix
            {
                rows_number: shape.1, columns_number: T::from(1),
                elements_indexes: unit_column_indexes, elements_values: unit_column_values
            });
            let unit_column = ExtendedMatrix { basic_matrix: basic_unit_column };
            let interim_inverse_column = l_matrix
                .naive_gauss_elimination(&unit_column).unwrap();
            let inverse_column = u_matrix
                .naive_gauss_elimination(&interim_inverse_column).unwrap();
            let all_inverse_column_values =
                inverse_column.basic_matrix.extract_all_elements_values();
            for i in 0..shape.0.into()
            {
                let current_inverse_column_element_value = extract_element_value(
                T::from(i), T::default(), &all_inverse_column_values);
                if current_inverse_column_element_value != V::default()
                {
                    inverse_matrix_indexes.push(T::from(i) * shape.1 + T::from(k));
                    inverse_matrix_values.push(current_inverse_column_element_value);
                }
            }
        }
        let basic_inverse_matrix = Box::new(NonSymmetricMatrix
            {
                rows_number: shape.0, columns_number: shape.1,
                elements_indexes: inverse_matrix_indexes, elements_values: inverse_matrix_values
            });
        let basic_inverse_matrix = basic_inverse_matrix.into_symmetric();
        Ok(ExtendedMatrix { basic_matrix: basic_inverse_matrix })
    }


    pub fn remove_zeros_rows_columns(&mut self) -> Vec<ZerosRowColumn<T>>
    {
        self.basic_matrix.remove_zeros_rows_columns()
    }


    pub fn remove_zeros_row(&mut self, row: T)
    {
        self.basic_matrix = self.basic_matrix.remove_zeros_row(row)
    }


    pub fn separate<'a>(self, positions: Vec<MatrixElementPosition<T>>)
        -> Result<SeparatedMatrix<T, V>, &'a str>
    {
        let shape = self.basic_matrix.get_shape();
        let all_elements_values =
            self.basic_matrix.extract_all_elements_values();
        let k_aa_rows_number = shape.0 - T::from(positions.len() as ElementsNumbers);
        let k_aa_columns_number = shape.1 - T::from(positions.len() as ElementsNumbers);
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
            ExtendedMatrix::create(k_aa_rows_number, k_aa_columns_number, k_aa_elements);
        let k_ab_rows_number = shape.0 - T::from(positions.len() as ElementsNumbers);
        let k_ab_columns_number = T::from(positions.len() as ElementsNumbers);
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
            ExtendedMatrix::create(k_ab_rows_number, k_ab_columns_number, k_ab_elements);
        let k_ba_rows_number = T::from(positions.len() as ElementsNumbers);
        let k_ba_columns_number = shape.1 - T::from(positions.len() as ElementsNumbers);
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
            ExtendedMatrix::create(k_ba_rows_number, k_ba_columns_number, k_ba_elements);
        let k_bb_rows_number = T::from(positions.len() as ElementsNumbers);
        let k_bb_columns_number = T::from(positions.len() as ElementsNumbers);
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
            ExtendedMatrix::create(k_bb_rows_number, k_bb_columns_number, k_bb_elements);
        Ok(SeparatedMatrix {
            k_aa: k_aa_matrix, k_ab: k_ab_matrix, k_ba: k_ba_matrix, k_bb: k_bb_matrix
        })
    }
}
