#[derive(Debug)]
pub enum StiffnessType
{
    Kuu,
    Kuth,
    Kthu,
    Kthth,
}


#[derive(Debug)]
pub struct StiffnessGroup<T>
{
    pub stiffness_type: StiffnessType,
    pub number_1: T,
    pub number_2: T,
    pub indexes: Vec<T>
}
