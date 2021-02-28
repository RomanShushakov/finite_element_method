use crate::fem::{StressStrainComponent, ForceComponent};
use crate::extended_matrix::ExtendedMatrix;

use crate::{ElementsNumbers, ElementsValues};


struct ElementIPStrains<V>
{
    strains_values: Vec<V>,
    strains_components: Vec<StressStrainComponent>,
}


struct ElementIPStresses<V>
{
    strains_values: Vec<V>,
    strains_components: Vec<StressStrainComponent>,
}


struct ElementIPForces<V>
{
    force_values: Vec<V>,
    force_component: Vec<ForceComponent>,
}


pub struct ElementAnalysisData<T, V>
{
    element_number: T,
    strains: Vec<ElementIPStrains<V>>,
    stresses: Vec<ElementIPStresses<V>>,
    forces: Vec<ElementIPForces<V>>,
}

//
// impl<T, V> ElementAnalysisData<T, V>
// {
//
// }


struct ElementsAnalysisResult<T, V>
{
    elements_analysis_data: Vec<ElementAnalysisData<T, V>>,
}
