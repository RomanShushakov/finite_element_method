use crate::fem::{DOFParameterData};
use crate::extended_matrix::ExtendedMatrix;

use crate::{ElementsNumbers, ElementsValues};

use std::ops::{Mul, Add, Sub, Div, Rem, SubAssign, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;


pub struct Displacements<T, V>
{
    pub displacements_values: Vec<V>,
    pub dof_parameters_data: Vec<DOFParameterData<T>>,
}


pub struct Reactions<T, V>
{
    pub reactions_values: Vec<V>,
    pub dof_parameters_data: Vec<DOFParameterData<T>>,
}


pub struct GlobalAnalysisResult<T, V>
{
    displacements: Displacements<T, V>,
    reactions: Reactions<T, V>,
}


impl<T, V> GlobalAnalysisResult<T, V>
    where T: Copy +PartialEq + Mul<Output = T> + Add<Output = T> + Sub<Output = T> +
             Div<Output = T> + Rem<Output = T> + Default + From<ElementsNumbers> +
             Into<ElementsNumbers> + Eq + Hash + SubAssign + Debug + PartialOrd + 'static,
          V: Copy + Default + Mul<Output = V> + Div<Output = V> + Sub<Output = V> +
             Add<Output = V> + From<ElementsValues> + Debug + PartialEq + AddAssign + MulAssign +
             SubAssign + Into<ElementsValues> + 'static,
{
    pub fn create(
        reactions_values: Vec<V>,
        reactions_dof_parameters_data: Vec<DOFParameterData<T>>,
        displacements_values: Vec<V>,
        displacements_dof_parameters_data: Vec<DOFParameterData<T>>,) -> Self
    {
        let reactions = Reactions { reactions_values,
            dof_parameters_data: reactions_dof_parameters_data };
        let displacements = Displacements { displacements_values,
            dof_parameters_data: displacements_dof_parameters_data };
        GlobalAnalysisResult { reactions, displacements }
    }


    pub fn extract_reactions(&self) -> Reactions<T, V>
    {
        let mut reactions_values = Vec::new();
        let mut dof_parameters_data = Vec::new();
        for (reaction_value, dof_parameter_data) in
            self.reactions.reactions_values.iter().zip(&self.reactions.dof_parameters_data)
        {
            reactions_values.push(*reaction_value);
            dof_parameters_data.push(*dof_parameter_data);
        }
        Reactions { reactions_values, dof_parameters_data }
    }


    pub fn extract_displacements(&self) -> Displacements<T, V>
    {
        let mut displacements_values = Vec::new();
        let mut dof_parameters_data = Vec::new();
        for (reaction_value, dof_parameter_data) in
            self.displacements.displacements_values.iter().zip(&self.displacements.dof_parameters_data)
        {
            displacements_values.push(*reaction_value);
            dof_parameters_data.push(*dof_parameter_data);
        }
        Displacements { displacements_values, dof_parameters_data }
    }
}
