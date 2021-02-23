use crate::fem::GLOBAL_DOF;

use std::slice::Iter;
use self::GlobalForceDisplacementComponent::*;


#[derive(PartialEq, Debug, Copy, Clone)]
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
    pub node_number: T,
    component: GlobalForceDisplacementComponent,
    value: V,
}


#[derive(Debug)]
pub struct ForceBC<T, V>
{
    pub number: T,
    pub force: Force<T, V>
}


impl<T, V> ForceBC<T, V>
{
    pub fn create(number: T, node_number: T, component: GlobalForceDisplacementComponent,
        value: V) -> Self
    {
        ForceBC { number, force: Force { node_number, component, value } }
    }


    pub fn update(&mut self, node_number: T, component: GlobalForceDisplacementComponent, value: V)
    {
        self.force.node_number = node_number;
        self.force.component = component;
        self.force.value = value;
    }
}


#[derive(Debug)]
pub struct Displacement<T, V>
{
    pub node_number: T,
    pub component: GlobalForceDisplacementComponent,
    pub value: V,
}


#[derive(Debug)]
pub struct DisplacementBC<T, V>
{
    pub number: T,
    pub displacement: Displacement<T, V>,
}


impl<T, V> DisplacementBC<T, V>
{
    pub fn create(number: T, node_number: T,
        component: GlobalForceDisplacementComponent, value: V) -> Self
    {
        DisplacementBC { number, displacement: Displacement { node_number, component, value } }
    }


    pub fn update(&mut self, node_number: T, component: GlobalForceDisplacementComponent, value: V)
    {
        self.displacement.node_number = node_number;
        self.displacement.component = component;
        self.displacement.value = value;
    }
}