use crate::ElementsNumbers;

use std::slice::Iter;
use self::StressStrainComponent::*;


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
