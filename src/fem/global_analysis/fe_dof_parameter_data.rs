use std::slice::Iter;
use std::ops::AddAssign;

use self::GlobalDOFParameter::*;


pub const GLOBAL_DOF: usize = 6;


pub fn global_dof<T>() -> T
    where T: AddAssign + From<u8>
{
    let mut global_dof = T::from(0u8);
    (0..GLOBAL_DOF).for_each(|_| global_dof += T::from(1u8));
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
    node_number: T,
    dof_parameter: GlobalDOFParameter,
}


impl<T> DOFParameterData<T>
    where T: Copy + PartialEq
{
    pub fn create(node_number: T, dof_parameter: GlobalDOFParameter) -> Self
    {
        DOFParameterData { node_number, dof_parameter }
    }


    pub fn node_number(&self) -> T
    {
        self.node_number
    }


    pub fn dof_parameter(&self) -> GlobalDOFParameter
    {
        self.dof_parameter
    }


    pub fn is_node_number_same(&self, node_number: T) -> bool
    {
        self.node_number == node_number
    }


    pub fn is_same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool
    {
        self.dof_parameter == dof_parameter && self.is_node_number_same(node_number)
    }


    pub fn update(&mut self, node_number: T, dof_parameter: GlobalDOFParameter)
    {
        self.node_number = node_number;
        self.dof_parameter = dof_parameter;
    }
}
