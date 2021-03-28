use crate::ElementsNumbers;

use std::slice::Iter;
use self::StressStrainComponent::*;
use crate::fem::EARComponentTrait;
use std::any::Any;


pub const STRESS_STRAIN_COMPONENTS_NUMBER: ElementsNumbers = 9;


#[derive(PartialEq, Debug, Copy, Clone)]
pub enum StressStrainComponent
{
    XX, XY, XZ,
    YX, YY, YZ,
    ZX, ZY, ZZ,
}


impl StressStrainComponent
{
    pub fn iterator() -> Iter<'static, StressStrainComponent>
     {
        const COMPONENTS: [StressStrainComponent; STRESS_STRAIN_COMPONENTS_NUMBER as usize] =
            [
                XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ,
            ];
        COMPONENTS.iter()
    }
}


impl StressStrainComponent
{
    pub fn as_str(&self) -> String
    {
        match self
        {
            StressStrainComponent::XX => String::from("XX"),
            StressStrainComponent::XY => String::from("XY"),
            StressStrainComponent::XZ => String::from("XZ"),
            StressStrainComponent::YX => String::from("YX"),
            StressStrainComponent::YY => String::from("YY"),
            StressStrainComponent::YZ => String::from("YZ"),
            StressStrainComponent::ZX => String::from("ZX"),
            StressStrainComponent::ZY => String::from("ZY"),
            StressStrainComponent::ZZ => String::from("ZZ"),
        }
    }
}


impl EARComponentTrait for StressStrainComponent
{
    fn as_any(&self) -> &dyn Any
    {
        self
    }


    fn same(&self, other: &Box<dyn EARComponentTrait>) -> bool
    {
        other
            .as_any()
            .downcast_ref::<StressStrainComponent>()
            .map_or(false, |component| self == component)
    }
}
