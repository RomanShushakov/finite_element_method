use std::collections::HashMap;

use extended_matrix::{SquareMatrix, FloatTrait, Vector};

use crate::fem::structs::{Props, Node, Truss};


pub const NODE_DOF: usize = 6;


pub struct FEM<V>
{
    props: Props<V>,
    stiffness_matrix: SquareMatrix<V>,
    displacements_vector: Vector<V>,
    forces_vector: Vector<V>,
    index: usize,
    indexes: Vec<usize>,
    imposed_constraints: Vec<bool>,
    imposed_forces: Vec<bool>,
    nodes: HashMap<u32, Node<V>>,
    truss_elements: HashMap<u32, Truss<V>>,
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
        let displacements_vector = Vector::create(&vec![V::from(0f32); nodes_number as usize * NODE_DOF]);
        let forces_vector = displacements_vector.clone();
        let index = 0;
        let indexes = (0..nodes_number as usize * NODE_DOF).collect::<Vec<usize>>();
        let imposed_constraints = vec![false; nodes_number as usize * NODE_DOF];
        let imposed_forces = vec![false; nodes_number as usize * NODE_DOF];

        let nodes = HashMap::new();
        let truss_elements = HashMap::new();

        FEM 
        { 
            props, stiffness_matrix, displacements_vector, forces_vector, index, indexes, imposed_constraints, 
            imposed_forces, nodes, truss_elements,
        }
    }


    pub(crate) fn get_stiffness_matrix(&self) -> &SquareMatrix<V>
    {
        &self.stiffness_matrix
    }


    pub(crate) fn get_displacements_vector(&self) -> &Vector<V>
    {
        &self.displacements_vector
    }


    pub(crate) fn get_forces_vector(&self) -> &Vector<V>
    {
        &self.forces_vector
    }


    pub(crate) fn get_index(&self) -> &usize
    {
        &self.index
    }


    pub(crate) fn get_mut_index(&mut self) -> &mut usize
    {
        &mut self.index
    }


    pub(crate) fn get_indexes(&self) -> &Vec<usize>
    {
        &self.indexes
    }


    pub(crate) fn get_imposed_constraints(&self) -> &Vec<bool>
    {
        &self.imposed_constraints
    }


    pub(crate) fn get_imposed_forces(&self) -> &Vec<bool>
    {
        &self.imposed_forces
    }


    pub(crate) fn get_nodes(&self) -> &HashMap<u32, Node<V>>
    {
        &self.nodes
    }


    pub(crate) fn get_mut_nodes(&mut self) -> &mut HashMap<u32, Node<V>>
    {
        &mut self.nodes
    }
}
