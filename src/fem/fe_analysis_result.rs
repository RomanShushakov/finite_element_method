use crate::fem::DOFParameterData;
use crate::extended_matrix::ExtendedMatrix;


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