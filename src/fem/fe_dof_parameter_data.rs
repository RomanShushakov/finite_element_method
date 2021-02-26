use crate::ElementsNumbers;

use std::slice::Iter;
use self::GlobalDOFParameter::*;


pub const GLOBAL_DOF: ElementsNumbers = 6;


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
        const PARAMETERS: [GlobalDOFParameter; GLOBAL_DOF as usize] =
            [
                X, Y, Z, ThX, ThY, ThZ,
            ];
        PARAMETERS.iter()
    }
}


#[derive(Debug, Copy, Clone)]
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
