use crate::fem::DOFParameterData;
use crate::extended_matrix::ExtendedMatrix;

use crate::{ElementsNumbers, ElementsValues};

use std::ops::{Mul, Add, Sub, Div, Rem, SubAssign, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;


struct Displacements<T, V>
{
    displacements_values: ExtendedMatrix<T, V>,
    dof_parameters_data: Vec<DOFParameterData<T>>,
}


struct Reactions<T, V>
{
    reactions_values: ExtendedMatrix<T, V>,
    dof_parameters_data: Vec<DOFParameterData<T>>,
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
        reactions_values: ExtendedMatrix<T, V>,
        reactions_dof_parameters_data: Vec<DOFParameterData<T>>,
        displacements_values: ExtendedMatrix<T, V>,
        displacements_dof_parameters_data: Vec<DOFParameterData<T>>,) -> Self
    {
        let reactions = Reactions { reactions_values,
            dof_parameters_data: reactions_dof_parameters_data };
        let displacements = Displacements { displacements_values,
            dof_parameters_data: displacements_dof_parameters_data };
        GlobalAnalysisResult { reactions, displacements }
    }


    pub fn show_reactions(&self)
    {
        println!("Reactions values:");
        self.reactions.reactions_values.show_matrix();
        println!("Reactions dof parameters data:");
        for parameter_data in &self.reactions.dof_parameters_data
        {
            println!("{:?}", parameter_data);
        }
    }


    pub fn show_displacements(&self)
    {
        println!("Displacements values:");
        self.displacements.displacements_values.show_matrix();
        println!("Displacements dof parameters data:");
        for parameter_data in &self.displacements.dof_parameters_data
        {
            println!("{:?}", parameter_data);
        }
    }
}
