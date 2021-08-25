use std::any::Any;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::hash::Hash;

use crate::fem::element_analysis::fe_stress_strain_components::StressStrainComponent;
use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::finite_elements::finite_element::FEType;


pub(crate) struct ElementStrains<V>
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


pub(crate) struct ElementStresses<V>
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


pub(crate) struct ElementForces<V>
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


pub struct NodalForces<V>
{
    forces_values: Vec<V>,
    forces_components: Vec<ForceComponent>,
}


impl<V> NodalForces<V>
{
    pub(crate) fn create(forces_values: Vec<V>, forces_components: Vec<ForceComponent>) -> Self
    {
        NodalForces { forces_values, forces_components }
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


pub struct ElementAnalysisData<T, V>
{
    strains: Option<ElementStrains<V>>,
    stresses: Option<ElementStresses<V>>,
    forces: Option<ElementForces<V>>,
    nodal_forces: Option<HashMap<T, NodalForces<V>>>,
}


impl<T, V> ElementAnalysisData<T, V>
    where T: Copy + PartialEq,
          V: Copy + PartialEq,
{
    pub(crate) fn create(strains: Option<ElementStrains<V>>, stresses: Option<ElementStresses<V>>,
        forces: Option<ElementForces<V>>, nodal_forces: Option<HashMap<T, NodalForces<V>>>) -> Self
    {
        ElementAnalysisData { strains, stresses, forces, nodal_forces }
    }


    pub(crate) fn extract_and_drop(self) -> (Option<ElementStrains<V>>, Option<ElementStresses<V>>,
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


    pub fn nodal_forces(&self) -> &Option<HashMap<T, NodalForces<V>>>
    {
        &self.nodal_forces
    }


    pub fn strains_values(&self) -> Option<&[V]>
    {
        if let Some(element_strains) = &self.strains
        {
            Some(element_strains.strains_values())
        }
        else
        {
            None
        }
    }


    pub fn strains_components(&self) -> Option<&[StressStrainComponent]>
    {
        if let Some(element_strains) = &self.strains
        {
            Some(element_strains.strains_components())
        }
        else
        {
            None
        }
    }

    pub fn stresses_values(&self) -> Option<&[V]>
    {
        if let Some(element_stresses) = &self.stresses
        {
            Some(element_stresses.stresses_values())
        }
        else
        {
            None
        }
    }


    pub fn stresses_components(&self) -> Option<&[StressStrainComponent]>
    {
        if let Some(element_stresses) = &self.stresses
        {
            Some(element_stresses.stresses_components())
        }
        else
        {
            None
        }
    }
}


pub struct ElementsAnalysisResult<T, V>
{
    elements_by_types: HashMap<FEType, Vec<T>>,
    elements_analysis_data: HashMap<T, ElementAnalysisData<T, V>>,
}


impl<T, V> ElementsAnalysisResult<T, V>
    where T: Debug + Eq + Hash
{
    pub(crate) fn create() -> Self
    {

        ElementsAnalysisResult { elements_by_types: HashMap::new(),
            elements_analysis_data: HashMap::new() }
    }


    pub(crate) fn add_to_types(&mut self, fe_type: FEType, number: T) -> Result<(), String>
    {
        if let Some(elements_numbers) = self.elements_by_types.get_mut(&fe_type)
        {
            if elements_numbers.iter().position(|element_number| *element_number == number)
                .is_some()
            {
                return Err(format!("ElementsAnalysisResult: Elements set {} does already contain \
                    number {:?}!", fe_type.as_str(), number));
            }
            elements_numbers.push(number);
        }
        else
        {
            self.elements_by_types.insert(fe_type, vec![number]);
        }

        Ok(())
    }


    pub(crate) fn add_to_analysis_data(&mut self, element_number: T,
        element_analysis_data: ElementAnalysisData<T, V>)
    {
        self.elements_analysis_data.insert(element_number, element_analysis_data);
    }


    pub fn elements_analysis_data(&self) -> &HashMap<T, ElementAnalysisData<T, V>>
    {
        &self.elements_analysis_data
    }


    pub fn elements_by_types(&self) -> &HashMap<FEType, Vec<T>>
    {
        &self.elements_by_types
    }
}


pub(crate) enum EARType
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


pub(super) trait EARComponentTrait: Any
{
    fn as_any(&self) -> &dyn Any;
    fn is_same(&self, other: &Box<dyn EARComponentTrait>) -> bool;
}
