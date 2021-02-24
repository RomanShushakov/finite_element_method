use crate::fem::GLOBAL_DOF;

use std::slice::Iter;
use self::GlobalDOFParameter::*;


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
        static PARAMETERS: [GlobalDOFParameter; GLOBAL_DOF as usize] =
            [
                X, Y, Z, ThX, ThY, ThZ,
            ];
        PARAMETERS.iter()
    }
}


#[derive(Debug)]
pub struct DOFParameterData<T>
{
    pub node_number: T,
    pub dof_parameter: GlobalDOFParameter,
}


impl<T> DOFParameterData<T>
    where T: PartialEq
{
    fn node_number_same(&self, node_number: T) -> bool
    {
        self.node_number == node_number
    }


    pub fn same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool
    {
        self.dof_parameter == dof_parameter && self.node_number_same(node_number)
    }
}


#[derive(Debug)]
pub struct Force<T, V>
{
    number: T,
    dof_parameter_data: DOFParameterData<T>,
    pub value: V,
}


impl<T, V> Force<T, V>
    where T: PartialEq
{
    pub fn create(number: T, node_number: T, dof_parameter: GlobalDOFParameter, value: V) -> Self
    {
        Force { number, dof_parameter_data: DOFParameterData { node_number, dof_parameter }, value }
    }


    pub fn update(&mut self, node_number: T, dof_parameter: GlobalDOFParameter, value: V)
    {
        self.dof_parameter_data.node_number = node_number;
        self.dof_parameter_data.dof_parameter = dof_parameter;
        self.value = value;
    }


    pub fn number_same(&self, number: T) -> bool
    {
        self.number == number
    }


    pub fn node_number_same(&self, node_number: T) -> bool
    {
        self.dof_parameter_data.node_number_same(node_number)
    }


    pub fn dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool
    {
        self.dof_parameter_data.same(dof_parameter, node_number)
    }
}


#[derive(Debug)]
pub struct Displacement<T, V>
{
    number: T,
    dof_parameter_data: DOFParameterData<T>,
    pub value: V,
}


impl<T, V> Displacement<T, V>
    where T: PartialEq
{
    pub fn create(number: T, node_number: T, dof_parameter: GlobalDOFParameter, value: V) -> Self
    {
        Displacement { number,
            dof_parameter_data: DOFParameterData { node_number, dof_parameter }, value }
    }


    pub fn update(&mut self, node_number: T, component: GlobalDOFParameter, value: V)
    {
        self.dof_parameter_data.node_number = node_number;
        self.dof_parameter_data.dof_parameter = component;
        self.value = value;
    }


    pub fn number_same(&self, number: T) -> bool
    {
        self.number == number
    }


    pub fn node_number_same(&self, node_number: T) -> bool
    {
        self.dof_parameter_data.node_number_same(node_number)
    }


    pub fn dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool
    {
        self.dof_parameter_data.same(dof_parameter, node_number)
    }
}