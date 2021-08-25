use std::any::Any;

use crate::fem::element_analysis::fe_element_analysis_result::EARComponentTrait;


#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum ForceComponent
{
    ForceX,
    ForceY,
    ForceZ,
    MomentX,
    MomentY,
    MomentZ,
}


impl ForceComponent
{
    pub fn as_str(&self) -> &'static str
    {
        match self
        {
            ForceComponent::ForceX => "ForceX",
            ForceComponent::ForceY => "ForceY",
            ForceComponent::ForceZ => "ForceZ",
            ForceComponent::MomentX => "MomentX",
            ForceComponent::MomentY => "MomentY",
            ForceComponent::MomentZ => "MomentZ",
        }
    }
}


impl EARComponentTrait for ForceComponent
{
    fn as_any(&self) -> &dyn Any
    {
        self
    }


    fn is_same(&self, other: &Box<dyn EARComponentTrait>) -> bool
    {
        other
            .as_any()
            .downcast_ref::<ForceComponent>()
            .map_or(false, |component| self == component)
    }
}
