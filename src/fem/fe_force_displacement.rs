use crate::fem::GLOBAL_DOF;

use std::slice::Iter;
use self::GlobalForceDisplacementComponent::*;


#[derive(PartialEq, Debug)]
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


#[derive(Debug)]
pub struct Force<T, V>
{
    pub force_number: T,
    pub node_number: T,
    pub component: GlobalForceDisplacementComponent,
    pub value: V,
}


impl<T, V> Force<T, V>
{
    pub fn create(force_number: T, node_number: T, component: GlobalForceDisplacementComponent,
        value: V) -> Self
    {
        Force { force_number, node_number, component, value }
    }


    pub fn update(&mut self, node_number: T,
        component: GlobalForceDisplacementComponent, value: V)
    {
        self.node_number = node_number;
        self.component = component;
        self.value = value;
    }
}


#[derive(Debug)]
pub struct Displacement<T, V>
{
    pub displacement_number: T,
    pub node_number: T,
    pub component: GlobalForceDisplacementComponent,
    pub value: V,
}


impl<T, V> Displacement<T, V>
{
    pub fn create(displacement_number: T, node_number: T,
        component: GlobalForceDisplacementComponent, value: V) -> Self
    {
        Displacement { displacement_number, node_number, component, value }
    }


    pub fn update(&mut self, node_number: T,
        component: GlobalForceDisplacementComponent, value: V)
    {
        self.node_number = node_number;
        self.component = component;
        self.value = value;
    }
}