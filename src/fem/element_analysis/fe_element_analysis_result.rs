use crate::fem::{StressStrainComponent, ForceComponent};
use crate::extended_matrix::ExtendedMatrix;

use crate::{ElementsNumbers, ElementsValues};


struct ElementStrain<V>
{
    strains_values: Vec<V>,
    strains_components: Vec<StressStrainComponent>,
}


struct ElementStress<V>
{
    strains_values: Vec<V>,
    strains_components: Vec<StressStrainComponent>,
}


struct ElementForce<V>
{
    force_value: V,
    force_component: ForceComponent,
}
