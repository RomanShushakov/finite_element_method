use std::slice::Iter;

use self::GlobalDOFParameter::*;


pub(crate) const GLOBAL_DOF: usize = 6;


pub(crate) fn global_dof() -> usize
{
    GLOBAL_DOF
}


#[derive(PartialEq, Debug, Copy, Clone, Hash, Eq)]
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
    pub(crate) fn iterator() -> Iter<'static, GlobalDOFParameter>
     {
        const PARAMETERS: [GlobalDOFParameter; GLOBAL_DOF] =
            [
                X, Y, Z, ThX, ThY, ThZ,
            ];
        PARAMETERS.iter()
    }
}


#[derive(Debug, Copy, Clone, PartialEq)]
pub struct DOFParameterData
{
    node_number: u32,
    dof_parameter: GlobalDOFParameter,
}


impl DOFParameterData
{
    pub(crate) fn create(node_number: u32, dof_parameter: GlobalDOFParameter) -> Self
    {
        DOFParameterData { node_number, dof_parameter }
    }


    pub fn copy_node_number(&self) -> u32
    {
        self.node_number
    }


    pub fn copy_dof_parameter(&self) -> GlobalDOFParameter
    {
        self.dof_parameter
    }


    pub(crate) fn is_node_number_same(&self, node_number: u32) -> bool
    {
        self.node_number == node_number
    }


    pub(crate) fn is_same(&self, dof_parameter: GlobalDOFParameter, node_number: u32) -> bool
    {
        self.dof_parameter == dof_parameter && self.is_node_number_same(node_number)
    }


    pub(crate) fn update(&mut self, node_number: u32, dof_parameter: GlobalDOFParameter)
    {
        self.node_number = node_number;
        self.dof_parameter = dof_parameter;
    }
}
