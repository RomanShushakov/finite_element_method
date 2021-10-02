use std::slice::Iter;
use std::ops::AddAssign;

use extended_matrix::matrix_element_position::MatrixElementPosition;

use self::StiffnessType::*;


const STIFFNESS_TYPES_NUMBER: usize = 4;


pub fn stiffness_types_number<T>() -> T
    where T: AddAssign + From<u8>
{
    let mut n = T::from(0u8);
    (0..STIFFNESS_TYPES_NUMBER).for_each(|_| n += T::from(1u8));
    n
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


pub struct StiffnessGroup<T>
{
    pub stiffness_type: StiffnessType,
    pub number_1: T,
    pub number_2: T,
    pub positions: Vec<MatrixElementPosition<T>>
}


#[derive(Hash, Eq, PartialEq)]
pub struct StiffnessGroupKey<T>
{
    pub stiffness_type: StiffnessType,
    pub number_1: T,
    pub number_2: T,
}
