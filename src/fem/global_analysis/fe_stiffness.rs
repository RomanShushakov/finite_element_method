use std::slice::Iter;

use extended_matrix::Position;

use self::StiffnessType::*;


const STIFFNESS_TYPES_NUMBER: usize = 4;


pub fn stiffness_types_number() -> usize
{
    STIFFNESS_TYPES_NUMBER
}


#[derive(Debug, Copy, Clone, PartialEq, Hash, Eq)]
pub enum StiffnessType
{
    Kuu,
    Kuth,
    Kthu,
    Kthth,
}


impl StiffnessType
{
    pub fn iterator() -> Iter<'static, StiffnessType>
     {
        const TYPES: [StiffnessType; STIFFNESS_TYPES_NUMBER] =
            [
                Kuu, Kuth, Kthu, Kthth,
            ];
        TYPES.iter()
    }
}


pub struct StiffnessGroup
{
    pub stiffness_type: StiffnessType,
    pub number_1: u32,
    pub number_2: u32,
    pub positions: Vec<Position>
}


#[derive(Hash, Eq, PartialEq)]
pub struct StiffnessGroupKey
{
    pub stiffness_type: StiffnessType,
    pub number_1: u32,
    pub number_2: u32,
}
