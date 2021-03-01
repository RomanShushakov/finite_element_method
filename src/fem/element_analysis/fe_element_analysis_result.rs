use crate::fem::{StressStrainComponent, ForceComponent};
use crate::extended_matrix::ExtendedMatrix;

use crate::{ElementsNumbers, ElementsValues};


#[derive(Debug)]
struct ElementStrains<V>
{
    strains_values: Vec<V>,
    strains_components: Vec<StressStrainComponent>,
}


#[derive(Debug)]
struct ElementStresses<V>
{
    stresses_values: Vec<V>,
    stresses_components: Vec<StressStrainComponent>,
}

#[derive(Debug)]
struct ElementForces<V>
{
    forces_values: Vec<V>,
    forces_components: Vec<ForceComponent>,
}


#[derive(Debug)]
pub struct ElementAnalysisData<T, V>
{
    element_number: T,
    strains: ElementStrains<V>,
    stresses: ElementStresses<V>,
    forces: ElementForces<V>,
}


impl<T, V> ElementAnalysisData<T, V>
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
}


#[derive(Debug)]
pub struct ElementsAnalysisResult<T, V>
{
    elements_analysis_data: Vec<ElementAnalysisData<T, V>>,
}


impl<T, V> ElementsAnalysisResult<T, V>
{
    pub fn create(elements_analysis_data: Vec<ElementAnalysisData<T, V>>) -> Self
    {
        ElementsAnalysisResult { elements_analysis_data }
    }
}
