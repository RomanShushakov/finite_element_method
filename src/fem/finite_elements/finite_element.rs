use std::hash::Hash;
use std::fmt::Debug;
use std::slice::Iter;
use std::collections::HashMap;
use std::any::Any;


use extended_matrix::{Matrix, FloatTrait};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::truss::truss2n1ip::Truss2n1ip;
use crate::fem::finite_elements::truss::truss2n2ip::Truss2n2ip;
use crate::fem::finite_elements::beam::beam2n1ipt::Beam2n1ipT;
use crate::fem::finite_elements::membrane::mem4n4ip::Mem4n4ip;
use crate::fem::finite_elements::plate::plate4n4ip::Plate4n4ip;
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
    Plate4n4ip,
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
            FEType::Mem4n4ip => "Mem4n4ip",
            FEType::Plate4n4ip => "Plate4n4ip",
        }
    }


    pub(super) fn iterator() -> Iter<'static, FEType>
    {
        const TYPES: [FEType; 5] =
            [
                Truss2n1ip, Truss2n2ip, Beam2n1ipT, Mem4n4ip, Plate4n4ip
            ];
        TYPES.iter()
    }
}


pub(crate) trait FiniteElementTrait<V>: Any
{
    fn update(
        &mut self, 
        nodes_numbers: Vec<u32>, 
        properties: Vec<V>,  
        nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    ) -> Result<(), String>;

    fn extract_stiffness_matrix(&self) -> Result<Matrix<V>, String>;

    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup>;

    fn is_node_belongs_to_element(&self, node_number: u32) -> bool;

    fn refresh(&mut self, nodes: &HashMap<u32, FENode<V>>, rel_tol: V, abs_tol: V) -> Result<(), String>;

    fn is_nodes_numbers_same(&self, nodes_numbers: Vec<u32>) -> bool;

    fn extract_element_analysis_data(
        &self, 
        global_displacements: &Displacements<V>,
        nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
    ) -> Result<ElementAnalysisData<V>, String>;

    fn copy_nodes_numbers(&self) -> Vec<u32>;

    fn extract_unique_elements_of_rotation_matrix(&self) -> Result<Vec<V>, String>;

    fn copy_properties(&self) -> Vec<V>;

    fn as_any(&self) -> &dyn Any;
}


struct FECreator<V>(V);


impl<V> FECreator<V>
    where V: FloatTrait<Output = V>
{
    fn create(
        fe_type: FEType, 
        nodes_numbers: Vec<u32>, 
        properties: Vec<V>, 
        ref_nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    ) 
        -> Result<Box<dyn FiniteElementTrait<V>>, String>
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
                            properties[0], 
                            properties[1],
                            Some(properties[2]), 
                            ref_nodes,
                            rel_tol,
                            abs_tol,
                        )?;

                        Ok(Box::new(truss_element))
                    }
                    else
                    {
                        let truss_element = Truss2n1ip::create(
                            nodes_numbers[0],
                            nodes_numbers[1],
                            properties[0], 
                            properties[1],
                            None, 
                            ref_nodes,
                            rel_tol,
                            abs_tol,
                        )?;

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
                            properties[0], 
                            properties[1],
                            Some(properties[2]), 
                            ref_nodes,
                            rel_tol,
                            abs_tol,
                        )?;

                        Ok(Box::new(truss_element))
                    }
                    else
                    {
                        let truss_element = Truss2n2ip::create(
                            nodes_numbers[0],
                            nodes_numbers[1],
                            properties[0], 
                            properties[1],
                            None, 
                            ref_nodes,
                            rel_tol,
                            abs_tol,
                        )?;

                        Ok(Box::new(truss_element))
                    }
                },
            FEType::Beam2n1ipT =>
                {
                    let beam_element = Beam2n1ipT::create(
                        nodes_numbers[0],
                        nodes_numbers[1],
                        properties[0], 
                        properties[1],
                        properties[2], 
                        properties[3],
                        properties[4], 
                        properties[5],
                        properties[6],
                        properties[7],
                        [properties[8], properties[9], properties[10]],
                        ref_nodes,
                        rel_tol, 
                        abs_tol,
                    )?;

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
                        ref_nodes,
                        rel_tol,
                        abs_tol,
                    )?;
                    
                    Ok(Box::new(membrane_element))
                },
            FEType::Plate4n4ip => 
                {
                    let plate_element = Plate4n4ip::create(
                        nodes_numbers[0],
                        nodes_numbers[1], 
                        nodes_numbers[2], 
                        nodes_numbers[3], 
                        properties[0], 
                        properties[1], 
                        properties[2],
                        properties[3], 
                        ref_nodes,
                        rel_tol,
                        abs_tol,
                    )?;
                    
                    Ok(Box::new(plate_element))
                }
        }
    }
}


pub(crate) struct FiniteElement<V>
{
    element_type: FEType,
    pub(crate) element: Box<dyn FiniteElementTrait<V>>,
}


impl<V> FiniteElement<V>
    where V: FloatTrait<Output = V>
{
    pub fn create(
        fe_type: FEType, 
        nodes_numbers: Vec<u32>, 
        properties: Vec<V>, 
        nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    ) 
        -> Result<Self, String>
    {
        let element = FECreator::create(
            fe_type, nodes_numbers, properties, nodes, rel_tol, abs_tol,
        )?;
        Ok(FiniteElement { element_type: fe_type, element })
    }


    pub fn update(
        &mut self, 
        nodes_numbers: Vec<u32>, 
        properties: Vec<V>, 
        nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    )
        -> Result<(), String>
    {
        self.element.update(nodes_numbers, properties, nodes, rel_tol, abs_tol)?;
        Ok(())
    }


    pub fn extract_stiffness_matrix(&self) -> Result<Matrix<V>, String>
    {
        let stiffness_matrix = self.element.extract_stiffness_matrix()?;
        Ok(stiffness_matrix)
    }


    pub fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup>
    {
        self.element.extract_stiffness_groups()
    }


    pub fn is_node_belong_element(&self, node_number: u32) -> bool
    {
        self.element.is_node_belongs_to_element(node_number)
    }


    pub fn refresh(&mut self, nodes: &HashMap<u32, FENode<V>>, rel_tol: V, abs_tol: V) -> Result<(), String>
    {
        self.element.refresh(nodes, rel_tol, abs_tol)?;
        Ok(())
    }


    pub fn is_type_same(&self, element_type: &FEType) -> bool
    {
        self.element_type == *element_type
    }


    pub fn is_nodes_numbers_same(&self, nodes_numbers: Vec<u32>) -> bool
    {
        self.element.is_nodes_numbers_same(nodes_numbers)
    }


    pub fn extract_element_analysis_data(
        &self, 
        global_displacements: &Displacements<V>,
        nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    )
        -> Result<ElementAnalysisData<V>, String>
    {
        let element_analysis_data = self.element
            .extract_element_analysis_data(global_displacements, nodes, rel_tol,
        )?;
        Ok(element_analysis_data)
    }


    pub fn copy_fe_type(&self) -> FEType
    {
        self.element_type
    }


    pub fn copy_nodes_numbers(&self) -> Vec<u32>
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
pub struct DeletedFEData<V>
{
    element_number: u32,
    element_type: FEType,
    nodes_numbers: Vec<u32>,
    properties: Vec<V>
}


impl<V> DeletedFEData<V>
    where V: FloatTrait<Output = V>
{
    pub(crate) fn create(element_number: u32, deleted_element: FiniteElement<V>) -> Self
    {
        let element_type = deleted_element.copy_fe_type();
        let nodes_numbers = deleted_element.copy_nodes_numbers();
        let properties = deleted_element.copy_properties();
        DeletedFEData { element_number, element_type, nodes_numbers, properties }
    }


    pub fn copy_number(&self) -> u32
    {
        self.element_number
    }


    pub fn copy_fe_type(&self) -> FEType
    {
        self.element_type
    }


    pub fn ref_nodes_numbers(&self) -> &[u32]
    {
        self.nodes_numbers.as_slice()
    }


    pub fn ref_properties(&self) -> &[V]
    {
        self.properties.as_slice()
    }
}
