use std::collections::HashMap;

use extended_matrix::{SquareMatrix, FloatTrait, Vector};

use crate::fem::structs::{Props, Node, NODE_DOF, Truss};


pub struct FEM<V>
{
    props: Props<V>,
    stiffness_matrix: SquareMatrix<V>,
    displacements_vector: Vector<V>,
    forces_vector: Vector<V>,
    nodes_count: usize,
    indexes: Vec<usize>,
    imposed_constraints: Vec<bool>,
    skyline: Vec<usize>,
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
        let nodes_count = 0;
        let indexes = (0..nodes_number as usize * NODE_DOF).collect::<Vec<usize>>();
        let imposed_constraints = vec![false; nodes_number as usize * NODE_DOF];
        let skyline = vec![0; nodes_number as usize * NODE_DOF];

        let nodes = HashMap::new();
        let truss_elements = HashMap::new();

        FEM 
        { 
            props, stiffness_matrix, displacements_vector, forces_vector, nodes_count, indexes, imposed_constraints,
            skyline, nodes, truss_elements,
        }
    }


    pub(crate) fn get_props(&self) -> &Props<V>
    {
        &self.props
    }


    pub(crate) fn get_stiffness_matrix(&self) -> &SquareMatrix<V>
    {
        &self.stiffness_matrix
    }


    pub(crate) fn get_mut_stiffness_matrix(&mut self) -> &mut SquareMatrix<V>
    {
        &mut self.stiffness_matrix
    }


    pub(crate) fn get_displacements_vector(&self) -> &Vector<V>
    {
        &self.displacements_vector
    }


    pub(crate) fn get_mut_displacements_vector(&mut self) -> &mut Vector<V>
    {
        &mut self.displacements_vector
    }


    pub(crate) fn get_forces_vector(&self) -> &Vector<V>
    {
        &self.forces_vector
    }


    pub(crate) fn get_mut_forces_vector(&mut self) -> &mut Vector<V>
    {
        &mut self.forces_vector
    }


    pub(crate) fn get_nodes_count(&self) -> &usize
    {
        &self.nodes_count
    }


    pub(crate) fn get_mut_nodes_count(&mut self) -> &mut usize
    {
        &mut self.nodes_count
    }


    pub(crate) fn get_indexes(&self) -> &Vec<usize>
    {
        &self.indexes
    }


    pub(crate) fn get_imposed_constraints(&self) -> &Vec<bool>
    {
        &self.imposed_constraints
    }


    pub(crate) fn get_mut_imposed_constraints(&mut self) -> &mut Vec<bool>
    {
        &mut self.imposed_constraints
    }


    pub(crate) fn get_skyline(&self) -> &Vec<usize>
    {
        &self.skyline
    }


    pub(crate) fn get_mut_skyline(&mut self) -> &mut Vec<usize>
    {
        &mut self.skyline
    }


    pub(crate) fn get_nodes(&self) -> &HashMap<u32, Node<V>>
    {
        &self.nodes
    }


    pub(crate) fn get_mut_nodes(&mut self) -> &mut HashMap<u32, Node<V>>
    {
        &mut self.nodes
    }


    pub(crate) fn get_truss_elements(&self) -> &HashMap<u32, Truss<V>>
    {
        &self.truss_elements
    }


    pub(crate) fn get_mut_truss_elements(&mut self) -> &mut HashMap<u32, Truss<V>>
    {
        &mut self.truss_elements
    }
}
