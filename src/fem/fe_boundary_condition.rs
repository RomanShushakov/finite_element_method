use crate::fem::{GlobalDOFParameter, DOFParameterData};
use std::fmt::Debug;


#[derive(Debug)]
struct Force<T, V>
{
    number: T,
    dof_parameter_data: DOFParameterData<T>,
    value: V,
}


impl<T, V> Force<T, V>
    where T: PartialEq
{
    fn create(number: T, node_number: T, dof_parameter: GlobalDOFParameter, value: V) -> Self
    {
        Force { number, dof_parameter_data: DOFParameterData { node_number, dof_parameter }, value }
    }
}


#[derive(Debug)]
struct Displacement<T, V>
{
    number: T,
    dof_parameter_data: DOFParameterData<T>,
    value: V,
}


impl<T, V> Displacement<T, V>
    where T: PartialEq
{
    fn create(number: T, node_number: T, dof_parameter: GlobalDOFParameter, value: V) -> Self
    {
        Displacement { number,
            dof_parameter_data: DOFParameterData { node_number, dof_parameter }, value }
    }
}


trait BCTrait<T, V>
{
    fn update(&mut self, node_number: T, dof_parameter: GlobalDOFParameter, value: V);
    fn number_same(&self, number: T) -> bool;
    fn node_number_same(&self, node_number: T) -> bool;
    fn dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool;
    fn extract_node_number(&self) -> T;
    fn extract_value(&self) -> V;
    fn extract_dof_parameter(&self) -> GlobalDOFParameter;
}


impl<T, V> BCTrait<T, V> for Force<T, V>
    where T: Copy + PartialEq,
          V: Copy,
{
    fn update(&mut self, node_number: T, dof_parameter: GlobalDOFParameter, value: V)
    {
        self.dof_parameter_data.node_number = node_number;
        self.dof_parameter_data.dof_parameter = dof_parameter;
        self.value = value;
    }


   fn number_same(&self, number: T) -> bool
    {
        self.number == number
    }


    fn node_number_same(&self, node_number: T) -> bool
    {
        self.dof_parameter_data.node_number_same(node_number)
    }


    fn dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool
    {
        self.dof_parameter_data.same(dof_parameter, node_number)
    }


    fn extract_node_number(&self) -> T
    {
        self.dof_parameter_data.node_number
    }


    fn extract_value(&self) -> V
    {
        self.value
    }


    fn extract_dof_parameter(&self) -> GlobalDOFParameter
    {
        self.dof_parameter_data.dof_parameter
    }
}


impl<T, V> BCTrait<T, V> for Displacement<T, V>
    where T: Copy + PartialEq,
          V: Copy,
{
    fn update(&mut self, node_number: T, dof_parameter: GlobalDOFParameter, value: V)
    {
        self.dof_parameter_data.node_number = node_number;
        self.dof_parameter_data.dof_parameter = dof_parameter;
        self.value = value;
    }


    fn number_same(&self, number: T) -> bool
    {
        self.number == number
    }


    fn node_number_same(&self, node_number: T) -> bool
    {
        self.dof_parameter_data.node_number_same(node_number)
    }


    fn dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool
    {
        self.dof_parameter_data.same(dof_parameter, node_number)
    }


    fn extract_node_number(&self) -> T
    {
        self.dof_parameter_data.node_number
    }


    fn extract_value(&self) -> V
    {
        self.value
    }


    fn extract_dof_parameter(&self) -> GlobalDOFParameter
    {
        self.dof_parameter_data.dof_parameter
    }
}


#[derive(Clone, Copy, PartialEq)]
pub enum BCType
{
    Force,
    Displacement,
}


impl BCType
{
    pub fn as_str(&self) -> String
    {
        match self
        {
            BCType::Force => String::from("Force"),
            BCType::Displacement => String::from("Displacement"),
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


pub struct BoundaryCondition<T, V>
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

    pub fn type_same(&self, bc_type: BCType) -> bool
    {
        self.bc_type == bc_type
    }


    pub fn number_same(&self, number: T) -> bool
    {
        self.boundary_condition.number_same(number)
    }


    pub fn node_number_same(&self, node_number: T) -> bool
    {
        self.boundary_condition.node_number_same(node_number)
    }


    pub fn dof_parameter_data_same(&self, dof_parameter: GlobalDOFParameter, node_number: T) -> bool
    {
        self.boundary_condition.dof_parameter_data_same(dof_parameter, node_number)
    }


    pub fn show(&self)
    {
        println!("{}, dof_param: {:?}, node: {:?}, value: {:?}", self.bc_type.as_str(),
                 self.boundary_condition.extract_dof_parameter(),
                 self.boundary_condition.extract_node_number(),
                 self.boundary_condition.extract_value())
    }
}
