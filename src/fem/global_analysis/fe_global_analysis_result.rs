use std::ops::{Mul, Add, Sub, Div, Rem, SubAssign, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;

use crate::fem::global_analysis::fe_dof_parameter_data::DOFParameterData;


pub struct Displacements<T, V>
{
    displacements_values: Vec<V>,
    dof_parameters_data: Vec<DOFParameterData<T>>,
}


impl<T, V> Displacements<T, V>
{
    fn create(displacements_values: Vec<V>, dof_parameters_data: Vec<DOFParameterData<T>>) -> Self
    {
        Displacements { displacements_values, dof_parameters_data }
    }


    pub fn displacements_values(&self) -> &[V]
    {
        self.displacements_values.as_slice()
    }


    pub fn dof_parameters_data(&self) -> &[DOFParameterData<T>]
    {
        self.dof_parameters_data.as_slice()
    }
}


pub struct Reactions<T, V>
{
    reactions_values: Vec<V>,
    dof_parameters_data: Vec<DOFParameterData<T>>,
}


impl<T, V> Reactions<T, V>
{
    fn create(reactions_values: Vec<V>, dof_parameters_data: Vec<DOFParameterData<T>>) -> Self
    {
        Reactions { reactions_values, dof_parameters_data }
    }

    pub fn reactions_values(&self) -> &[V]
    {
        self.reactions_values.as_slice()
    }


    pub fn dof_parameters_data(&self) -> &[DOFParameterData<T>]
    {
        self.dof_parameters_data.as_slice()
    }
}


pub struct GlobalAnalysisResult<T, V>
{
    displacements: Displacements<T, V>,
    reactions: Reactions<T, V>,
}


impl<T, V> GlobalAnalysisResult<T, V>
    where T: Copy +PartialEq + Mul<Output = T> + Add<Output = T> + Sub<Output = T> +
             Div<Output = T> + Rem<Output = T> + Eq + Hash + SubAssign + Debug +
             PartialOrd + 'static,
          V: Copy + Mul<Output = V> + Div<Output = V> + Sub<Output = V> +
             Add<Output = V> + Debug + PartialEq + AddAssign + MulAssign +
             SubAssign + Into<f64> + 'static,
{
    pub fn create(
        reactions_values: Vec<V>,
        reactions_dof_parameters_data: Vec<DOFParameterData<T>>,
        displacements_values: Vec<V>,
        displacements_dof_parameters_data: Vec<DOFParameterData<T>>,) -> Self
    {
        let reactions = Reactions::create(reactions_values,
            reactions_dof_parameters_data);
        let displacements = Displacements::create(displacements_values,
            displacements_dof_parameters_data);
        GlobalAnalysisResult { reactions, displacements }
    }


    pub fn reactions(&self) -> &Reactions<T, V>
    {
        &self.reactions
    }


    pub fn displacements(&self) -> &Displacements<T, V>
    {
        &self.displacements
    }
}
