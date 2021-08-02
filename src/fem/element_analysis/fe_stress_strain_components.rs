use std::slice::Iter;
use std::any::Any;

use crate::fem::element_analysis::fe_element_analysis_result::EARComponentTrait;

use self::StressStrainComponent::*;


pub const STRESS_STRAIN_COMPONENTS_NUMBER: usize = 9;


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
        const COMPONENTS: [StressStrainComponent; STRESS_STRAIN_COMPONENTS_NUMBER] =
            [
                XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ,
            ];
        COMPONENTS.iter()
    }
}


impl StressStrainComponent
{
    pub fn as_str(&self) -> &'static str
    {
        match self
        {
            StressStrainComponent::XX => "XX",
            StressStrainComponent::XY => "XY",
            StressStrainComponent::XZ => "XZ",
            StressStrainComponent::YX => "YX",
            StressStrainComponent::YY => "YY",
            StressStrainComponent::YZ => "YZ",
            StressStrainComponent::ZX => "ZX",
            StressStrainComponent::ZY => "ZY",
            StressStrainComponent::ZZ => "ZZ",
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
