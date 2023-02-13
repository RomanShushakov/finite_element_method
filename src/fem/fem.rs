use std::collections::HashMap;

use extended_matrix::{SquareMatrix, FloatTrait, Vector};

use crate::fem::structs::{Props, Node, NODE_DOF, Truss, Beam, Plate};


pub struct FEM<V>
{
    props: Props<V>,
    stiffness_matrix: SquareMatrix<V>,
    displacements_vector: Vector<V>,
    forces_vector: Vector<V>,
    nodes_count: usize,
    indexes: Vec<usize>,
    index_node_number_map: HashMap<usize, u32>,
    imposed_constraints: Vec<bool>,
    nodes: HashMap<u32, Node<V>>,
    truss_elements: HashMap<u32, Truss<V>>,
    beam_elements: HashMap<u32, Beam<V>>,
    plate_elements: HashMap<u32, Plate<V>>,
}


impl<V> FEM<V>
    where V: FloatTrait<Output = V>
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
        let index_node_number_map = HashMap::new();
        let imposed_constraints = vec![false; nodes_number as usize * NODE_DOF];

        let nodes = HashMap::new();
        let truss_elements = HashMap::new();
        let beam_elements = HashMap::new();
        let plate_elements = HashMap::new();

        FEM 
        { 
            props, stiffness_matrix, displacements_vector, forces_vector, nodes_count, indexes, index_node_number_map,
            imposed_constraints, nodes, truss_elements, beam_elements, plate_elements,
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


    pub(crate) fn get_index_node_number_map(&self) -> &HashMap<usize, u32>
    {
        &self.index_node_number_map
    }


    pub(crate) fn get_mut_index_node_number_map(&mut self) -> &mut HashMap<usize, u32>
    {
        &mut self.index_node_number_map
    }


    pub(crate) fn get_imposed_constraints(&self) -> &Vec<bool>
    {
        &self.imposed_constraints
    }


    pub(crate) fn get_mut_imposed_constraints(&mut self) -> &mut Vec<bool>
    {
        &mut self.imposed_constraints
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


    pub(crate) fn get_beam_elements(&self) -> &HashMap<u32, Beam<V>>
    {
        &self.beam_elements
    }


    pub(crate) fn get_mut_beam_elements(&mut self) -> &mut HashMap<u32, Beam<V>>
    {
        &mut self.beam_elements
    }


    pub(crate) fn get_plate_elements(&self) -> &HashMap<u32, Plate<V>>
    {
        &self.plate_elements
    }


    pub(crate) fn get_mut_plate_elements(&mut self) -> &mut HashMap<u32, Plate<V>>
    {
        &mut self.plate_elements
    }


    pub fn reset(&mut self, nodes_number: u32)
    {
        self.stiffness_matrix = SquareMatrix::create(
            nodes_number as usize * NODE_DOF, &Vec::new(),
        );
        self.displacements_vector = Vector::create(&vec![V::from(0f32); nodes_number as usize * NODE_DOF]);
        self.forces_vector = self.displacements_vector.clone();
        self.nodes_count = 0;
        self.indexes = (0..nodes_number as usize * NODE_DOF).collect::<Vec<usize>>();
        self.index_node_number_map = HashMap::new();
        self.imposed_constraints = vec![false; nodes_number as usize * NODE_DOF];

        self.nodes = HashMap::new();
        self.truss_elements = HashMap::new();
        self.beam_elements = HashMap::new();
        self.plate_elements = HashMap::new();
    }


    pub fn get_truss_rotation_matrix_elements(&self, number: u32) -> Result<[V; 9], String>
    {
        self.check_truss_element_exist(number)?;

        let truss = self.get_truss_elements()
            .get(&number)
            .ok_or(format!("Truss element {number} is absent!"))?;

        Ok(truss.get_rotation_matrix_elements())
    }


    pub fn get_beam_rotation_matrix_elements(&self, number: u32) -> Result<[V; 9], String>
    {
        self.check_beam_element_exist(number)?;

        let beam = self.get_beam_elements()
            .get(&number)
            .ok_or(format!("Beam element {number} is absent!"))?;

        Ok(beam.get_rotation_matrix_elements())
    }


    pub fn get_plate_rotation_matrix_elements(&self, number: u32) -> Result<[V; 9], String>
    {
        self.check_plate_element_exist(number)?;

        let plate = self.get_plate_elements()
            .get(&number)
            .ok_or(format!("Plate element {number} is absent!"))?;

        Ok(plate.get_rotation_matrix_elements())
    }
}
