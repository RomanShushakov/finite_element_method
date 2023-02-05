use std::any::Any;
use std::collections::HashMap;
use std::fmt::Debug;

use crate::fem::element_analysis::fe_stress_strain_components::StressStrainComponent;
use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::finite_elements::finite_element::FEType;


#[derive(Debug, Clone)]
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


    pub fn ref_strains_values(&self) -> &[V]
    {
        self.strains_values.as_slice()
    }


    pub fn ref_strains_components(&self) -> &[StressStrainComponent]
    {
        self.strains_components.as_slice()
    }
}


#[derive(Debug, Clone)]
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


    pub fn ref_stresses_values(&self) -> &[V]
    {
        self.stresses_values.as_slice()
    }


    pub fn ref_stresses_components(&self) -> &[StressStrainComponent]
    {
        self.stresses_components.as_slice()
    }
}


#[derive(Debug, Clone)]
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


    pub fn ref_forces_values(&self) -> &[V]
    {
        self.forces_values.as_slice()
    }


    pub fn ref_forces_components(&self) -> &[ForceComponent]
    {
        self.forces_components.as_slice()
    }
}


#[derive(Debug, Clone)]
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


    pub fn ref_forces_values(&self) -> &[V]
    {
        self.forces_values.as_slice()
    }


    pub fn ref_forces_components(&self) -> &[ForceComponent]
    {
        self.forces_components.as_slice()
    }
}


#[derive(Debug, Clone)]
pub struct ElementAnalysisData<V>
{
    optional_strains: Option<ElementStrains<V>>,
    optional_stresses: Option<ElementStresses<V>>,
    optional_forces: Option<ElementForces<V>>,
    optional_nodal_forces: Option<HashMap<u32, NodalForces<V>>>,
}


impl<V> ElementAnalysisData<V>
    where V: Copy + PartialEq,
{
    pub(crate) fn create(
        optional_strains: Option<ElementStrains<V>>,
        optional_stresses: Option<ElementStresses<V>>,
        optional_forces: Option<ElementForces<V>>,
        optional_nodal_forces: Option<HashMap<u32, NodalForces<V>>>
    ) 
        -> Self
    {
        ElementAnalysisData { optional_strains, optional_stresses, optional_forces, optional_nodal_forces }
    }


    pub(crate) fn extract_and_drop(self) -> (Option<ElementStrains<V>>, Option<ElementStresses<V>>,
        Option<ElementForces<V>>)
    {
        (self.optional_strains, self.optional_stresses,
         self.optional_forces)
    }


    pub fn optional_ref_forces_values(&self) -> Option<&[V]>
    {
        if let Some(element_forces) = &self.optional_forces
        {
            Some(element_forces.ref_forces_values())
        }
        else
        {
            None
        }
    }


    pub fn optional_ref_forces_components(&self) -> Option<&[ForceComponent]>
    {
        if let Some(element_forces) = &self.optional_forces
        {
            Some(element_forces.ref_forces_components())
        }
        else
        {
            None
        }
    }


    pub fn ref_optional_nodal_forces(&self) -> &Option<HashMap<u32, NodalForces<V>>>
    {
        &self.optional_nodal_forces
    }


    pub fn optional_ref_strains_values(&self) -> Option<&[V]>
    {
        if let Some(element_strains) = &self.optional_strains
        {
            Some(element_strains.ref_strains_values())
        }
        else
        {
            None
        }
    }


    pub fn optional_ref_strains_components(&self) -> Option<&[StressStrainComponent]>
    {
        if let Some(element_strains) = &self.optional_strains
        {
            Some(element_strains.ref_strains_components())
        }
        else
        {
            None
        }
    }

    pub fn optional_ref_stresses_values(&self) -> Option<&[V]>
    {
        if let Some(element_stresses) = &self.optional_stresses
        {
            Some(element_stresses.ref_stresses_values())
        }
        else
        {
            None
        }
    }


    pub fn optional_ref_stresses_components(&self) -> Option<&[StressStrainComponent]>
    {
        if let Some(element_stresses) = &self.optional_stresses
        {
            Some(element_stresses.ref_stresses_components())
        }
        else
        {
            None
        }
    }
}


#[derive(Debug, Clone)]
pub struct ElementsAnalysisResult<V>
{
    elements_by_types: HashMap<FEType, Vec<u32>>,
    elements_analysis_data: HashMap<u32, ElementAnalysisData<V>>,
}


impl<V> ElementsAnalysisResult<V>
{
    pub(crate) fn create() -> Self
    {

        ElementsAnalysisResult { elements_by_types: HashMap::new(), elements_analysis_data: HashMap::new() }
    }


    pub(crate) fn add_to_types(&mut self, fe_type: FEType, number: u32) -> Result<(), String>
    {
        if let Some(elements_numbers) = self.elements_by_types.get_mut(&fe_type)
        {
            if elements_numbers.iter().position(|element_number| *element_number == number).is_some()
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


    pub(crate) fn add_to_analysis_data(&mut self, element_number: u32, element_analysis_data: ElementAnalysisData<V>)
    {
        self.elements_analysis_data.insert(element_number, element_analysis_data);
    }


    pub fn ref_elements_analysis_data(&self) -> &HashMap<u32, ElementAnalysisData<V>>
    {
        &self.elements_analysis_data
    }


    pub fn ref_elements_by_types(&self) -> &HashMap<FEType, Vec<u32>>
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
