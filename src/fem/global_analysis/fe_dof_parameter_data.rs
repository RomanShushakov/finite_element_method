use std::slice::Iter;
use std::ops::AddAssign;

use extended_matrix::one::One;

use self::GlobalDOFParameter::*;



pub const GLOBAL_DOF: usize = 6;


pub fn global_dof<T>() -> T
    where T: One + Default + AddAssign
{
    let mut global_dof = T::default();
    (0..GLOBAL_DOF).for_each(|_| global_dof += T::one());
    global_dof
}


#[derive(PartialEq, Debug, Copy, Clone)]
pub enum GlobalDOFParameter
{
    X,
    Y,
    Z,
    ThX,
    ThY,
    ThZ,
}


impl GlobalDOFParameter
{
    pub fn iterator() -> Iter<'static, GlobalDOFParameter>
     {
        const PARAMETERS: [GlobalDOFParameter; GLOBAL_DOF] =
            [
                X, Y, Z, ThX, ThY, ThZ,
            ];
        PARAMETERS.iter()
    }
}


#[derive(Debug, Copy, Clone, PartialEq)]
pub struct DOFParameterData<T>
{
    pub node_number: T,
    pub dof_parameter: GlobalDOFParameter,
}


impl<T> DOFParameterData<T>
    where T: PartialEq
{
    pub fn node_number_same(&self, node_number: T) -> bool
    {
        self.node_number == node_number
    }


    pub fn same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool
    {
        self.dof_parameter == dof_parameter && self.node_number_same(node_number)
    }
}
