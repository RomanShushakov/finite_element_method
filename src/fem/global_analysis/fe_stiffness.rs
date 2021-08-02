use std::slice::Iter;
use std::ops::AddAssign;

use extended_matrix::one::One;
use extended_matrix::basic_matrix::basic_matrix::MatrixElementPosition;

use crate::ElementsNumbers;

use self::StiffnessType::*;



pub const STIFFNESS_TYPES_NUMBER: ElementsNumbers = 4;


pub fn stiffness_types_number<T>() -> T
    where T: Default + One + AddAssign
{
    let mut n = T::default();
    (0..STIFFNESS_TYPES_NUMBER).for_each(|_| n += T::one());
    n
}


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
