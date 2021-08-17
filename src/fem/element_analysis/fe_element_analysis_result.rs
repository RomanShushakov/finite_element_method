use std::any::Any;

use crate::fem::element_analysis::fe_stress_strain_components::StressStrainComponent;
use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;


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
        self.strains.clone()
    }


    pub fn extract_stresses(&self) -> ElementStresses<V>
    {
        self.stresses.clone()
    }


    pub fn extract_forces(&self) -> ElementForces<V>
    {
        self.forces.clone()
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
        self.elements_analysis_data.clone()
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
    pub fn as_str(&self) -> &'static str
    {
        match self
        {
            EARType::Stress => "Stress",
            EARType::Strain => "Strain",
            EARType::Force => "Force",
        }
    }
}


pub trait EARComponentTrait: Any
{
    fn as_any(&self) -> &dyn Any;
    fn is_same(&self, other: &Box<dyn EARComponentTrait>) -> bool;
}
