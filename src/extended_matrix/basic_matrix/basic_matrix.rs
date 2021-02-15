use std::fmt::Debug;
use std::any::Any;
use std::collections::HashMap;


#[derive(Debug)]
pub enum BasicMatrixType
{
    Symmetric,
    NonSymmetric
}


pub struct Shape<T>(pub T, pub T);


#[derive(Debug, Hash, PartialEq, Eq)]
pub struct MatrixElementPosition<T>
{
    pub row: T,
    pub column: T
}


#[derive(Debug)]
pub struct ZerosRowColumn<T>
{
    pub row: T,
    pub column: T,
}


pub trait BasicMatrixClone<T, V>
{
    fn clone_box(&self) -> Box<dyn BasicMatrix<T, V>>;
}


impl<T, V, W> BasicMatrixClone<T, V> for W
    where W: BasicMatrix<T, V> + Clone + 'static,
{
    fn clone_box(&self) -> Box<dyn BasicMatrix<T, V>>
    {
        Box::new(self.clone())
    }
}


impl<T, V> Clone for Box<dyn BasicMatrix<T, V>>
{
    fn clone(&self) -> Box<dyn BasicMatrix<T, V>>
    {
        self.clone_box()
    }
}


pub trait BasicMatrix<T, V>: BasicMatrixClone<T, V>
{
    // fn create_element_value(&mut self, requested_index: T, new_value: V);
    fn read_element_value(&self, row: T, column: T) -> Result<V, &str>;
    // fn update_element_value(&mut self, row: T, column: T, new_value: V) -> Result<(), &str>;
    // fn delete_element_value(&mut self, row: T, column: T) -> Result<(), &str>;
    fn extract_all_elements_values(&self) -> HashMap<MatrixElementPosition<T>, V>;
    fn get_shape(&self) -> Shape<T>;
    fn transpose(&mut self);
    fn multiply_by_number(&mut self, number: V);
    fn into_symmetric(self) -> Box<dyn BasicMatrix<T, V>>;
    fn define_type(&self) -> BasicMatrixType;
    fn as_any(&self) -> &dyn Any;
    fn remove_zeros_rows_columns(&mut self) -> Vec<ZerosRowColumn<T>>;
}
