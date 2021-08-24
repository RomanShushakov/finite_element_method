use std::any::Any;
use std::collections::HashMap;

use crate::fem::element_analysis::fe_stress_strain_components::StressStrainComponent;
use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::finite_elements::finite_element::FEType;


#[derive(Debug, PartialEq)]
pub struct ElementStrains<V>
{
    strains_values: Vec<V>,
    strains_components: Vec<StressStrainComponent>,
}


impl<V> ElementStrains<V>
{
    pub fn create(strains_values: Vec<V>, strains_components: Vec<StressStrainComponent>) -> Self
    {
        ElementStrains { strains_values, strains_components }
    }


    pub fn strains_values(&self) -> &[V]
    {
        self.strains_values.as_slice()
    }


    pub fn strains_components(&self) -> &[StressStrainComponent]
    {
        self.strains_components.as_slice()
    }
}


#[derive(Debug, PartialEq)]
pub struct ElementStresses<V>
{
    stresses_values: Vec<V>,
    stresses_components: Vec<StressStrainComponent>,
}


impl<V> ElementStresses<V>
{
    pub fn create(stresses_values: Vec<V>, stresses_components: Vec<StressStrainComponent>) -> Self
    {
        ElementStresses { stresses_values, stresses_components }
    }


    pub fn stresses_values(&self) -> &[V]
    {
        self.stresses_values.as_slice()
    }


    pub fn stresses_components(&self) -> &[StressStrainComponent]
    {
        self.stresses_components.as_slice()
    }
}


#[derive(Debug, PartialEq)]
pub struct ElementForces<V>
{
    forces_values: Vec<V>,
    forces_components: Vec<ForceComponent>,
}


impl<V> ElementForces<V>
{
    pub fn create(forces_values: Vec<V>, forces_components: Vec<ForceComponent>) -> Self
    {
        ElementForces { forces_values, forces_components }
    }


    pub fn forces_values(&self) -> &[V]
    {
        self.forces_values.as_slice()
    }


    pub fn forces_components(&self) -> &[ForceComponent]
    {
        self.forces_components.as_slice()
    }
}


#[derive(Debug, PartialEq)]
pub struct ElementAnalysisData<V>
{
    strains: Option<ElementStrains<V>>,
    stresses: Option<ElementStresses<V>>,
    forces: Option<ElementForces<V>>,
}


impl<V> ElementAnalysisData<V>
    where V: Copy + PartialEq,
{
    pub fn create(strains: Option<ElementStrains<V>>, stresses: Option<ElementStresses<V>>,
        forces: Option<ElementForces<V>>) -> Self
    {
        ElementAnalysisData { strains, stresses, forces }
    }


    pub fn extract_and_drop(self) -> (Option<ElementStrains<V>>, Option<ElementStresses<V>>,
        Option<ElementForces<V>>)
    {
        (self.strains, self.stresses, self.forces)
    }


    pub fn forces_values(&self) -> Option<&[V]>
    {
        if let Some(element_forces) = &self.forces
        {
            Some(element_forces.forces_values())
        }
        else
        {
            None
        }
    }


    pub fn forces_components(&self) -> Option<&[ForceComponent]>
    {
        if let Some(element_forces) = &self.forces
        {
            Some(element_forces.forces_components())
        }
        else
        {
            None
        }
    }
}


pub struct ElementsAnalysisResult<T, V>
{
    elements_analysis_data: HashMap<T, ElementAnalysisData<V>>,     // Hashmap { element_number: ElementAnalysisData }
    analyzed_elements_types: HashMap<FEType, Vec<T>>,               // Hashmap { FEType: Vec<element_number, ...> }
}


impl<T, V> ElementsAnalysisResult<T, V>
{
    pub fn create(elements_analysis_data: HashMap<T, ElementAnalysisData<V>>,
        analyzed_elements_types: HashMap<FEType, Vec<T>>) -> Self
    {
        ElementsAnalysisResult { elements_analysis_data, analyzed_elements_types }
    }


    pub fn analyzed_elements_types(&self) -> &HashMap<FEType, Vec<T>>
    {
        &self.analyzed_elements_types
    }


    pub fn elements_analysis_data(&self) -> &HashMap<T, ElementAnalysisData<V>>
    {
        &self.elements_analysis_data
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
