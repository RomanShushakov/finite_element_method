use std::collections::HashMap;
use extended_matrix::{SquareMatrix, FloatTrait, Vector};

use crate::fem::structs::{Props, Node};


pub const NODE_DOF: usize = 6;

pub struct FEM<V>
{
    props: Props<V>,
    pub stiffness_matrix: SquareMatrix<V>,
    pub displacements: Vector<V>,
    pub forces: Vector<V>,
    nodes: HashMap<u32, Node<V>>,
}


impl<V> FEM<V>
    where V: FloatTrait
{
    pub fn create(rel_tol: V, abs_tol: V, nodes_number: u32) -> Self
    {
        let props = Props::create(rel_tol, abs_tol, nodes_number);
        let stiffness_matrix = SquareMatrix::create(
            nodes_number as usize * NODE_DOF, &Vec::new(),
        );
        let displacements = Vector::create(&vec![V::from(0f32); nodes_number as usize * NODE_DOF]);
        let forces = displacements.clone();
        let nodes = HashMap::new();

        FEM { props, stiffness_matrix, displacements, forces, nodes }
    }
}
