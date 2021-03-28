use crate::fem::{StressStrainComponent, ForceComponent};

use std::any::Any;


#[derive(Debug, Clone, PartialEq)]
pub struct ElementStrains<V>
{
    pub strains_values: Vec<V>,
    pub strains_components: Vec<StressStrainComponent>,
}


#[derive(Debug, Clone, PartialEq)]
pub struct ElementStresses<V>
{
    pub stresses_values: Vec<V>,
    pub stresses_components: Vec<StressStrainComponent>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ElementForces<V>
{
    pub forces_values: Vec<V>,
    pub forces_components: Vec<ForceComponent>,
}


#[derive(Debug, Clone, PartialEq)]
pub struct ElementAnalysisData<T, V>
{
    element_number: T,
    strains: ElementStrains<V>,
    stresses: ElementStresses<V>,
    forces: ElementForces<V>,
}


impl<T, V> ElementAnalysisData<T, V>
    where T: Copy + PartialEq,
          V: Copy + PartialEq,
{
    pub fn create(element_number: T,
        strains_values: Vec<V>, strains_components: Vec<StressStrainComponent>,
        stresses_values: Vec<V>, stresses_components: Vec<StressStrainComponent>,
        forces_values: Vec<V>, forces_components: Vec<ForceComponent>) -> Self
    {
        let strains = ElementStrains { strains_values, strains_components };
        let stresses = ElementStresses { stresses_values, stresses_components };
        let forces = ElementForces { forces_values, forces_components };
        ElementAnalysisData { element_number, strains, stresses, forces }
    }


    pub fn number_same(&self, number: T) -> bool
    {
        self.element_number == number
    }


    pub fn extract_element_number(&self) -> T
    {
        self.element_number
    }


    pub fn extract_strains(&self) -> ElementStrains<V>
    {
        self.strains.to_owned()
    }


    pub fn extract_stresses(&self) -> ElementStresses<V>
    {
        self.stresses.to_owned()
    }


    pub fn extract_forces(&self) -> ElementForces<V>
    {
        self.forces.to_owned()
    }
}


#[derive(Debug)]
pub struct ElementsAnalysisResult<T, V>
{
    elements_analysis_data: Vec<ElementAnalysisData<T, V>>,
}


impl<T, V> ElementsAnalysisResult<T, V>
    where T: Copy,
          V: Copy,
{
    pub fn create(elements_analysis_data: Vec<ElementAnalysisData<T, V>>) -> Self
    {
        ElementsAnalysisResult { elements_analysis_data }
    }


    pub fn extract_elements_analysis_data(&self) -> Vec<ElementAnalysisData<T, V>>
    {
        self.elements_analysis_data.to_owned()
    }
}


pub enum EARType
{
    Stress,
    Strain,
    Force,
}


impl EARType
{
    pub fn as_str(&self) -> String
    {
        match self
        {
            EARType::Stress => String::from("Stress"),
            EARType::Strain => String::from("Strain"),
            EARType::Force => String::from("Force"),
        }
    }
}


pub trait EARComponentTrait: Any
{
    fn as_any(&self) -> &dyn Any;
    fn same(&self, other: &Box<dyn EARComponentTrait>) -> bool;
}
