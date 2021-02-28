use crate::extended_matrix::MatrixElementPosition;
use crate::ElementsNumbers;

use std::slice::Iter;
use self::StiffnessType::*;


pub const STIFFNESS_TYPES_NUMBER: ElementsNumbers = 4;


#[derive(Debug, Copy, Clone, PartialEq)]
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
        const TYPES: [StiffnessType; STIFFNESS_TYPES_NUMBER as usize] =
            [
                Kuu, Kuth, Kthu, Kthth,
            ];
        TYPES.iter()
    }
}


pub struct StiffnessGroup<T>
{
    pub stiffness_type: StiffnessType,
    pub number_1: T,
    pub number_2: T,
    pub positions: Vec<MatrixElementPosition<T>>
}
