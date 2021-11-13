use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::slice::Iter;
use std::collections::HashMap;

use extended_matrix::extended_matrix::ExtendedMatrix;

use crate::my_float::MyFloatTrait;

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::truss::truss2n1ip::Truss2n1ip;
use crate::fem::finite_elements::truss::truss2n2ip::Truss2n2ip;
use crate::fem::finite_elements::beam::beam2n1ipt::Beam2n1ipT;
use crate::fem::finite_elements::membrane::mem4n4ip::Mem4n4ip;
use crate::fem::global_analysis::fe_stiffness::StiffnessGroup;
use crate::fem::element_analysis::fe_element_analysis_result::ElementAnalysisData;
use crate::fem::global_analysis::fe_global_analysis_result::Displacements;


use self::FEType::*;


#[derive(Copy, Clone, PartialEq, Debug, Eq, Hash)]
pub enum FEType
{
    Truss2n1ip,
    Truss2n2ip,
    Beam2n1ipT,
    Mem4n4ip,
}


impl FEType
{
    pub fn as_str(&self) -> &'static str
    {
        match self
        {
            FEType::Truss2n1ip => "Truss2n1ip",
            FEType::Truss2n2ip => "Truss2n2ip",
            FEType::Beam2n1ipT => "Beam2n1ipT",
            &FEType::Mem4n4ip => "Mem4n4ip",
        }
    }


    pub(super) fn iterator() -> Iter<'static, FEType>
    {
        const TYPES: [FEType; 4] =
            [
                Truss2n1ip, Truss2n2ip, Beam2n1ipT, Mem4n4ip,
            ];
        TYPES.iter()
    }
}


pub(super) trait FiniteElementTrait<T, V>
{
    fn update(&mut self, nodes_numbers: Vec<T>, properties: Vec<V>, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<(), String>;
    fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>;
    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>;
    fn is_node_belongs_to_element(&self, node_number: T) -> bool;
    fn refresh(&mut self, tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<(), String>;
    fn is_nodes_numbers_same(&self, nodes_numbers: Vec<T>) -> bool;
    fn extract_element_analysis_data(&self, global_displacements: &Displacements<T, V>,
        tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<ElementAnalysisData<T, V>, String>;
    fn copy_nodes_numbers(&self) -> Vec<T>;
    fn extract_unique_elements_of_rotation_matrix(&self) -> Result<Vec<V>, String>;
    fn copy_properties(&self) -> Vec<V>;
}


struct FECreator<T, V>(T, V);


impl<T, V> FECreator<T, V>
    where T: Copy + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash + SubAssign +
             Debug + Mul<Output = T> + PartialOrd + Add<Output = T> + AddAssign + From<u8> +
             Ord + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             PartialEq + Debug + AddAssign + MulAssign + SubAssign + MyFloatTrait + PartialOrd +
             Into<f64> + From<f32> + MyFloatTrait<Other = V> + 'static,
{
    fn create(fe_type: FEType, nodes_numbers: Vec<T>, properties: Vec<V>, tolerance: V,
        ref_nodes: &HashMap<T, FENode<V>>) -> Result<Box<dyn FiniteElementTrait<T, V>>, String>
    {
        match fe_type
        {
            FEType::Truss2n1ip =>
                {
                    if properties.len() == 3
                    {
                        let truss_element = Truss2n1ip::create(
                            nodes_numbers[0],
                            nodes_numbers[1],
                            properties[0], properties[1],
                            Some(properties[2]), tolerance, ref_nodes)?;

                        Ok(Box::new(truss_element))
                    }
                    else
                    {
                        let truss_element = Truss2n1ip::create(
                            nodes_numbers[0],
                            nodes_numbers[1],
                            properties[0], properties[1],
                            None, tolerance, ref_nodes)?;

                        Ok(Box::new(truss_element))
                    }
                },
            FEType::Truss2n2ip =>
                {
                    if properties.len() == 3
                    {
                        let truss_element = Truss2n2ip::create(
                            nodes_numbers[0],
                            nodes_numbers[1],
                            properties[0], properties[1],
                            Some(properties[2]), tolerance, ref_nodes)?;

                        Ok(Box::new(truss_element))
                    }
                    else
                    {
                        let truss_element = Truss2n2ip::create(
                            nodes_numbers[0],
                            nodes_numbers[1],
                            properties[0], properties[1],
                            None, tolerance, ref_nodes)?;

                        Ok(Box::new(truss_element))
                    }
                },
            FEType::Beam2n1ipT =>
                {
                    let beam_element = Beam2n1ipT::create(
                    nodes_numbers[0],
                    nodes_numbers[1],
                    properties[0], properties[1],
                    properties[2], properties[3],
                    properties[4], properties[5],
                    properties[6],
                    properties[7],
                    [properties[8], properties[9], properties[10]],
                    tolerance, ref_nodes)?;

                    Ok(Box::new(beam_element))
                },
            FEType::Mem4n4ip => 
                {
                    let membrane_element = Mem4n4ip::create(
                        nodes_numbers[0],
                        nodes_numbers[1], 
                        nodes_numbers[2], 
                        nodes_numbers[3], 
                        properties[0], 
                        properties[1], 
                        properties[2], 
                        tolerance, ref_nodes)?;
                    
                    Ok(Box::new(membrane_element))
                }
        }
    }
}


pub(crate) struct FiniteElement<T, V>
{
    element_type: FEType,
    element: Box<dyn FiniteElementTrait<T, V>>,
}


impl<T, V> FiniteElement<T, V>
    where T: Copy + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash + SubAssign +
             Debug + Mul<Output = T> + PartialOrd + Add<Output = T> + AddAssign + From<u8> +
             Ord + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             PartialEq + Debug + AddAssign + MulAssign + SubAssign + MyFloatTrait + PartialOrd +
             Into<f64> + From<f32> + MyFloatTrait<Other = V> + 'static,
{
    pub fn create(fe_type: FEType, nodes_numbers: Vec<T>, properties: Vec<V>, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<Self, String>
    {
        let element = FECreator::create(fe_type,
            nodes_numbers, properties, tolerance, nodes)?;
        Ok(FiniteElement { element_type: fe_type, element })
    }


    pub fn update(&mut self, nodes_numbers: Vec<T>, properties: Vec<V>, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        self.element.update(nodes_numbers, properties, tolerance, nodes)?;
        Ok(())
    }


    pub fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
    {
        let stiffness_matrix = self.element.extract_stiffness_matrix()?;
        Ok(stiffness_matrix)
    }


    pub fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>
    {
        self.element.extract_stiffness_groups()
    }


    pub fn is_node_belong_element(&self, node_number: T) -> bool
    {
        self.element.is_node_belongs_to_element(node_number)
    }


    pub fn refresh(&mut self, tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        self.element.refresh(tolerance, nodes)?;
        Ok(())
    }


    pub fn is_type_same(&self, element_type: &FEType) -> bool
    {
        self.element_type == *element_type
    }


    pub fn is_nodes_numbers_same(&self, nodes_numbers: Vec<T>) -> bool
    {
        self.element.is_nodes_numbers_same(nodes_numbers)
    }


    pub fn extract_element_analysis_data(&self, global_displacements: &Displacements<T, V>,
        tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<ElementAnalysisData<T, V>, String>
    {
        let element_analysis_data = self.element.extract_element_analysis_data(
            global_displacements, tolerance, nodes)?;
        Ok(element_analysis_data)
    }


    pub fn copy_fe_type(&self) -> FEType
    {
        self.element_type
    }


    pub fn copy_nodes_numbers(&self) -> Vec<T>
    {
        self.element.copy_nodes_numbers()
    }


    pub fn extract_unique_elements_of_rotation_matrix(&self) -> Result<Vec<V>, String>
    {
        self.element.extract_unique_elements_of_rotation_matrix()
    }


    pub fn copy_properties(&self) -> Vec<V>
    {
        self.element.copy_properties()
    }
}


#[derive(Debug, Clone)]
pub struct DeletedFEData<T, V>
{
    element_number: T,
    element_type: FEType,
    nodes_numbers: Vec<T>,
    properties: Vec<V>
}


impl<T, V> DeletedFEData<T, V>
    where T: Copy + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash + SubAssign +
             Debug + Mul<Output = T> + PartialOrd + Add<Output = T> + AddAssign + From<u8> +
             Ord + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             PartialEq + Debug + AddAssign + MulAssign + SubAssign + MyFloatTrait + PartialOrd +
             Into<f64> + From<f32> + MyFloatTrait<Other = V> + 'static,
{
    pub(crate) fn create(element_number: T, deleted_element: FiniteElement<T, V>) -> Self
    {
        let element_type = deleted_element.copy_fe_type();
        let nodes_numbers = deleted_element.copy_nodes_numbers();
        let properties = deleted_element.copy_properties();
        DeletedFEData { element_number, element_type, nodes_numbers, properties }
    }


    pub fn copy_number(&self) -> T
    {
        self.element_number
    }


    pub fn copy_fe_type(&self) -> FEType
    {
        self.element_type
    }


    pub fn ref_nodes_numbers(&self) -> &[T]
    {
        self.nodes_numbers.as_slice()
    }


    pub fn ref_properties(&self) -> &[V]
    {
        self.properties.as_slice()
    }
}
