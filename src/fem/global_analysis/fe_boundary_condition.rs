use std::fmt::Debug;

use crate::fem::global_analysis::fe_dof_parameter_data::{GlobalDOFParameter, DOFParameterData};


struct Force<T, V>
{
    number: T,
    dof_parameter_data: DOFParameterData<T>,
    value: V,
}


impl<T, V> Force<T, V>
    where T: Copy + PartialEq
{
    fn create(number: T, node_number: T, dof_parameter: GlobalDOFParameter, value: V) -> Self
    {
        let dof_parameter_data =
            DOFParameterData::create(node_number, dof_parameter);
        Force { number, dof_parameter_data, value }
    }
}


struct Displacement<T, V>
{
    number: T,
    dof_parameter_data: DOFParameterData<T>,
    value: V,
}


impl<T, V> Displacement<T, V>
    where T: Copy + PartialEq
{
    fn create(number: T, node_number: T, dof_parameter: GlobalDOFParameter, value: V) -> Self
    {
        let dof_parameter_data = DOFParameterData::create(
            node_number, dof_parameter);
        Displacement { number, dof_parameter_data, value }
    }
}


trait BCTrait<T, V>
{
    fn update(&mut self, node_number: T, dof_parameter: GlobalDOFParameter, value: V);
    fn is_number_same(&self, number: T) -> bool;
    fn is_node_number_same(&self, node_number: T) -> bool;
    fn is_dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool;
    fn copy_node_number(&self) -> T;
    fn copy_value(&self) -> V;
    fn copy_dof_parameter(&self) -> GlobalDOFParameter;
    fn copy_number(&self) -> T;
}


impl<T, V> BCTrait<T, V> for Force<T, V>
    where T: Copy + PartialEq,
          V: Copy,
{
    fn update(&mut self, node_number: T, dof_parameter: GlobalDOFParameter, value: V)
    {
        self.dof_parameter_data.update(node_number, dof_parameter);
        self.value = value;
    }


   fn is_number_same(&self, number: T) -> bool
    {
        self.number == number
    }


    fn is_node_number_same(&self, node_number: T) -> bool
    {
        self.dof_parameter_data.is_node_number_same(node_number)
    }


    fn is_dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool
    {
        self.dof_parameter_data.is_same(dof_parameter, node_number)
    }


    fn copy_node_number(&self) -> T
    {
        self.dof_parameter_data.copy_node_number()
    }


    fn copy_value(&self) -> V
    {
        self.value
    }


    fn copy_dof_parameter(&self) -> GlobalDOFParameter
    {
        self.dof_parameter_data.copy_dof_parameter()
    }


    fn copy_number(&self) -> T
    {
        self.number
    }
}


impl<T, V> BCTrait<T, V> for Displacement<T, V>
    where T: Copy + PartialEq,
          V: Copy,
{
    fn update(&mut self, node_number: T, dof_parameter: GlobalDOFParameter, value: V)
    {
        self.dof_parameter_data.update(node_number, dof_parameter);
        self.value = value;
    }


    fn is_number_same(&self, number: T) -> bool
    {
        self.number == number
    }


    fn is_node_number_same(&self, node_number: T) -> bool
    {
        self.dof_parameter_data.is_node_number_same(node_number)
    }


    fn is_dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T)
                                  -> bool
    {
        self.dof_parameter_data.is_same(dof_parameter, node_number)
    }


    fn copy_node_number(&self) -> T
    {
        self.dof_parameter_data.copy_node_number()
    }


    fn copy_value(&self) -> V
    {
        self.value
    }


    fn copy_dof_parameter(&self) -> GlobalDOFParameter
    {
        self.dof_parameter_data.copy_dof_parameter()
    }


    fn copy_number(&self) -> T
    {
        self.number
    }
}


#[derive(Copy, Clone, PartialEq, Debug)]
pub enum BCType
{
    Force,
    Displacement,
}


impl BCType
{
    pub fn as_str(&self) -> &'static str
    {
        match self
        {
            BCType::Force => "Force",
            BCType::Displacement => "Displacement",
        }
    }
}


struct BCCreator<T, V>(T, V);


impl<T, V> BCCreator<T, V>
    where T: PartialEq + Copy + 'static,
          V: Copy + 'static
{
    fn create(bc_type: BCType, number: T, node_number: T,
        dof_parameter: GlobalDOFParameter, value: V) -> Box<dyn BCTrait<T, V>>
    {
        match bc_type
        {
            BCType::Force => Box::new(Force::create(number, node_number, dof_parameter, value)),
            BCType::Displacement => Box::new(Displacement::create(
                number, node_number, dof_parameter, value)),
        }
    }
}


pub(crate) struct BoundaryCondition<T, V>
{
    bc_type: BCType,
    boundary_condition: Box<dyn BCTrait<T, V>>
}


impl<T, V> BoundaryCondition<T, V>
    where T: Copy + PartialEq + Debug + 'static,
          V: Copy + Debug + 'static
{
    pub fn create(bc_type: BCType, number: T, node_number: T,
        dof_parameter: GlobalDOFParameter, value: V) -> Self
    {
        let boundary_condition =
            BCCreator::create(bc_type.clone(), number, node_number, dof_parameter, value);
        BoundaryCondition { bc_type, boundary_condition }
    }


    pub fn update(&mut self, node_number: T, dof_parameter: GlobalDOFParameter, value: V)
    {
        self.boundary_condition.update(node_number, dof_parameter, value)
    }


    pub fn is_type_same(&self, bc_type: BCType) -> bool
    {
        self.bc_type == bc_type
    }


    pub fn is_number_same(&self, number: T) -> bool
    {
        self.boundary_condition.is_number_same(number)
    }


    pub fn is_node_number_same(&self, node_number: T) -> bool
    {
        self.boundary_condition.is_node_number_same(node_number)
    }


    pub fn is_dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T)
        -> bool
    {
        self.boundary_condition.is_dof_parameter_data_same(dof_parameter, node_number)
    }


    pub fn copy_bc_type(&self) -> BCType
    {
        self.bc_type
    }


    pub fn copy_value(&self) -> V
    {
        self.boundary_condition.copy_value()
    }


    pub fn copy_number(&self) -> T
    {
        self.boundary_condition.copy_number()
    }


    pub fn copy_node_number(&self) -> T
    {
        self.boundary_condition.copy_node_number()
    }


    pub fn copy_dof_parameter(&self) -> GlobalDOFParameter
    {
        self.boundary_condition.copy_dof_parameter()
    }
}


#[derive(Debug, Clone)]
pub struct DeletedBCData<T, V>
{
    bc_type: BCType,
    number: T,
    node_number: T,
    dof_parameter: GlobalDOFParameter,
    value: V
}


impl<T, V> DeletedBCData<T, V>
    where T: Copy + PartialEq + Debug + 'static,
          V: Copy + Debug + 'static,
{
    pub(crate) fn create(deleted_bc: BoundaryCondition<T, V>) -> Self
    {
        let bc_type = deleted_bc.copy_bc_type();
        let number = deleted_bc.copy_number();
        let node_number = deleted_bc.copy_node_number();
        let dof_parameter = deleted_bc.copy_dof_parameter();
        let value = deleted_bc.copy_value();
        DeletedBCData { bc_type, number, node_number, dof_parameter, value }
    }


    pub fn copy_bc_type(&self) -> BCType
    {
        self.bc_type
    }


    pub fn copy_number(&self) -> T
    {
        self.number
    }


    pub fn copy_node_number(&self) -> T
    {
        self.node_number
    }


    pub fn copy_dof_parameter(&self) -> GlobalDOFParameter
    {
        self.dof_parameter
    }


    pub fn copy_value(&self) -> V
    {
        self.value
    }
}
