use std::any::Any;
use std::collections::HashMap;
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
    fe_type: FEType,
    strains: Option<ElementStrains<V>>,
    stresses: Option<ElementStresses<V>>,
    forces: Option<ElementForces<V>>,
    nodal_forces: Option<HashMap<T, NodalForces<V>>>,
}


impl<T, V> ElementAnalysisData<T, V>
    where T: Copy + PartialEq,
          V: Copy + PartialEq,
{
    pub(crate) fn create(fe_type: FEType, strains: Option<ElementStrains<V>>,
        stresses: Option<ElementStresses<V>>, forces: Option<ElementForces<V>>,
        nodal_forces: Option<HashMap<T, NodalForces<V>>>) -> Self
    {
        ElementAnalysisData { fe_type, strains, stresses, forces, nodal_forces }
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


    pub fn fe_type(&self) -> &FEType
    {
        &self.fe_type
    }
}


pub struct ElementsAnalysisResult<T, V>
{
    elements_analysis_data: HashMap<T, ElementAnalysisData<T, V>>,
    elements_with_strains_components: Option<HashMap<StressStrainComponent, Vec<T>>>,
    elements_with_stresses_components: Option<HashMap<StressStrainComponent, Vec<T>>>,
    elements_with_forces_components: Option<HashMap<ForceComponent, Vec<T>>>,
}


impl<T, V> ElementsAnalysisResult<T, V>
    where T: Debug + Eq + Hash
{
    pub(crate) fn create() -> Self
    {

        ElementsAnalysisResult { elements_analysis_data: HashMap::new(),
            elements_with_strains_components: None, elements_with_stresses_components: None,
            elements_with_forces_components: None }
    }


    pub(crate) fn add_to_analysis_data(&mut self, element_number: T,
        element_analysis_data: ElementAnalysisData<T, V>)
    {
        self.elements_analysis_data.insert(element_number, element_analysis_data);
    }


    pub(crate) fn add_to_stresses_strains_components(&mut self, ear_type: EARType,
        stress_strain_component: StressStrainComponent, element_number: T) -> Result<(), String>
    {
        let elements_with_stresses_strains_components = match ear_type
            {
                EARType::Strain => Ok(self.elements_with_strains_components.as_mut()),
                EARType::Stress => Ok(self.elements_with_stresses_components.as_mut()),
                EARType::Force => Err("ElementsAnalysisResult: Incorrect element analysis result \
                    type!"),
            }?;

        if let Some(stresses_strains_components) =
        elements_with_stresses_strains_components
        {
            if let Some(element_numbers) = stresses_strains_components
                .get_mut(&stress_strain_component)
            {
                if element_numbers.iter().position(|number| *number == element_number)
                    .is_some()
                {
                    return Err(format!("ElementsAnalysisResult: {} component {} does already \
                        contain element with number {:?}", ear_type.as_str(),
                        stress_strain_component.as_str(), element_number));
                }
                else
                {
                    element_numbers.push(element_number);
                }
            }
            else
            {
                stresses_strains_components.insert(stress_strain_component,
                    vec![element_number]);
            }
        }
        else
        {
            let mut elements_with_stresses_strains_components =
                HashMap::new();
            elements_with_stresses_strains_components.insert(stress_strain_component,
                vec![element_number]);

            match ear_type
            {
                EARType::Strain =>
                    {
                        self.elements_with_strains_components =
                            Some(elements_with_stresses_strains_components)
                    },
                EARType::Stress =>
                    {
                        self.elements_with_stresses_components =
                            Some(elements_with_stresses_strains_components)
                    }
                _ => ()
            }
        }
        Ok(())
    }


    pub(crate) fn add_to_forces_components(&mut self, force_component: ForceComponent,
        element_number: T) -> Result<(), String>
    {
        if let Some(elements_with_forces_components) =
            self.elements_with_forces_components.as_mut()
        {
            if let Some(element_numbers) = elements_with_forces_components
                .get_mut(&force_component)
            {
                if element_numbers.iter().position(|number| *number == element_number)
                    .is_some()
                {
                    return Err(format!("ElementsAnalysisResult: Force component {} does already \
                        contain element with number {:?}", force_component.as_str(),
                        element_number));
                }
                else
                {
                    element_numbers.push(element_number);
                }
            }
            else
            {
                elements_with_forces_components.insert(force_component, vec![element_number]);
            }
        }
        else
        {
            let mut elements_with_forces_components = HashMap::new();
            elements_with_forces_components.insert(force_component, vec![element_number]);
            self.elements_with_forces_components = Some(elements_with_forces_components);
        }
        Ok(())
    }


    pub fn elements_analysis_data(&self) -> &HashMap<T, ElementAnalysisData<T, V>>
    {
        &self.elements_analysis_data
    }


    pub fn elements_with_strains_components(&self)
        -> &Option<HashMap<StressStrainComponent, Vec<T>>>
    {
        &self.elements_with_strains_components
    }


    pub fn elements_with_stresses_components(&self)
        -> &Option<HashMap<StressStrainComponent, Vec<T>>>
    {
        &self.elements_with_stresses_components
    }


    pub fn elements_with_forces_components(&self) -> &Option<HashMap<ForceComponent, Vec<T>>>
    {
        &self.elements_with_forces_components
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
