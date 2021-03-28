use crate::fem::EARComponentTrait;
use std::any::Any;

#[derive(Debug, Clone, PartialEq, Copy)]
pub enum ForceComponent
{
    Axial
}


impl ForceComponent
{
    pub fn as_str(&self) -> String
    {
        match self
        {
            ForceComponent::Axial => String::from("Axial"),
        }
    }
}


impl EARComponentTrait for ForceComponent
{
    fn as_any(&self) -> &dyn Any
    {
        self
    }


    fn same(&self, other: &Box<dyn EARComponentTrait>) -> bool
    {
        other
            .as_any()
            .downcast_ref::<ForceComponent>()
            .map_or(false, |component| self == component)
    }
}
