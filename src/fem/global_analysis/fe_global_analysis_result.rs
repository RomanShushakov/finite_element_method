use std::ops::{Mul, Add, Sub, Div, SubAssign, AddAssign, MulAssign};
use std::fmt::Debug;

use crate::fem::global_analysis::fe_dof_parameter_data::DOFParameterData;


#[derive(Debug, Clone)]
pub struct Displacements<V>
{
    displacements_values: Vec<V>,
    dof_parameters_data: Vec<DOFParameterData>,
}


impl<V> Displacements<V>
{
    fn create(displacements_values: Vec<V>, dof_parameters_data: Vec<DOFParameterData>) -> Self
    {
        Displacements { displacements_values, dof_parameters_data }
    }


    pub(crate) fn displacements_values(&self) -> &[V]
    {
        self.displacements_values.as_slice()
    }


    pub(crate) fn dof_parameters_data(&self) -> &[DOFParameterData]
    {
        self.dof_parameters_data.as_slice()
    }
}


#[derive(Debug, Clone)]
struct Reactions<V>
{
    reactions_values: Vec<V>,
    dof_parameters_data: Vec<DOFParameterData>,
}


impl<V> Reactions<V>
{
    fn create(reactions_values: Vec<V>, dof_parameters_data: Vec<DOFParameterData>) -> Self
    {
        Reactions { reactions_values, dof_parameters_data }
    }

    fn reactions_values(&self) -> &[V]
    {
        self.reactions_values.as_slice()
    }


    fn dof_parameters_data(&self) -> &[DOFParameterData]
    {
        self.dof_parameters_data.as_slice()
    }
}


#[derive(Debug, Clone)]
pub struct GlobalAnalysisResult<V>
{
    displacements: Displacements<V>,
    reactions: Reactions<V>,
}


impl<V> GlobalAnalysisResult<V>
    where V: Copy + Mul<Output = V> + Div<Output = V> + Sub<Output = V> + Add<Output = V> + Debug + PartialEq + 
             AddAssign + MulAssign + SubAssign + Into<f64> + 'static,
{
    pub(crate) fn create(
        reactions_values: Vec<V>,
        reactions_dof_parameters_data: Vec<DOFParameterData>,
        displacements_values: Vec<V>,
        displacements_dof_parameters_data: Vec<DOFParameterData>,) -> Self
    {
        let reactions = Reactions::create(reactions_values,
            reactions_dof_parameters_data);
        let displacements = Displacements::create(displacements_values,
            displacements_dof_parameters_data);
        GlobalAnalysisResult { reactions, displacements }
    }


    pub fn reactions_values(&self) -> &[V]
    {
        self.reactions.reactions_values()
    }


    pub fn reactions_dof_parameters_data(&self) -> &[DOFParameterData]
    {
        self.reactions.dof_parameters_data()
    }


    pub fn displacements(&self) -> &Displacements<V>
    {
        &self.displacements
    }


    pub fn displacements_values(&self) -> &[V]
    {
        self.displacements.displacements_values()
    }


    pub fn displacements_dof_parameters_data(&self) -> &[DOFParameterData]
    {
        self.displacements.dof_parameters_data()
    }
}
