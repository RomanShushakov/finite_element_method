use crate::{FeNode, GlobalCoordinates};


pub struct Truss2n2ip<'a, T, V>
{
    pub number: T,
    pub node_1: &'a FeNode<T, V>,
    pub node_2: &'a FeNode<T, V>,
    pub young_modulus: V,
    pub area: V,
    pub area_2: Option<V>,
}