use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use std::collections::HashMap;


#[derive(Debug, Clone, PartialEq)]
pub struct NodalForces<V>
{
    forces_values: Vec<V>,
    forces_components: Vec<ForceComponent>,
}


impl<V> NodalForces<V>
{
    pub fn create(forces_values: Vec<V>, forces_components: Vec<ForceComponent>) -> Self
    {
        NodalForces { forces_values, forces_components }
    }


    pub fn forces_values(&self) -> &[V]
    {
        self.forces_values.as_slice()
    }


    pub fn forces_components(&self) -> &[ForceComponent]
    {
        self.forces_components.as_slice()
    }
}


pub struct BeamElementNodalForces<T, V>
{
    beam_element_nodal_forces: HashMap<T, NodalForces<V>>
}


impl<T, V> BeamElementNodalForces<T, V>
{
    pub fn create(beam_element_nodal_forces: HashMap<T, NodalForces<V>>) -> Self
    {
        BeamElementNodalForces { beam_element_nodal_forces }
    }


    pub fn beam_element_nodal_forces(&self) -> &HashMap<T, NodalForces<V>>
    {
        &self.beam_element_nodal_forces
    }
}
