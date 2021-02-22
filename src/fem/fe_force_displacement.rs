use crate::fem::GLOBAL_DOF;

use std::slice::Iter;
use self::GlobalForceDisplacementComponent::*;


pub enum GlobalForceDisplacementComponent
{
    X,
    Y,
    Z,
    ThX,
    ThY,
    ThZ,
}

impl GlobalForceDisplacementComponent
{
    pub fn iterator() -> Iter<'static, GlobalForceDisplacementComponent>
     {
        static COMPONENTS: [GlobalForceDisplacementComponent; GLOBAL_DOF as usize] =
            [
                X, Y, Z, ThX, ThY, ThZ,
            ];
        COMPONENTS.iter()
    }
}


pub struct Force<T, V>
{
    pub node_number: T,
    pub component: GlobalForceDisplacementComponent,
    pub value: V,
}


pub struct Displacement<T, V>
{
    pub node_number: T,
    pub component: GlobalForceDisplacementComponent,
    pub value: V,
}
