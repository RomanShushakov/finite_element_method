// use std::fmt::Debug;

use extended_matrix::{FloatTrait, BasicOperationsTrait, Position};

use crate::fem::FEM;
use crate::fem::structs::{Truss, TRUSS_NODE_DOF, NODE_DOF};


enum TrussError
{
    Number(u32),
    SameNodes(u32, u32),
}


impl TrussError
{
    fn compose_error_message(&self) -> String
    {
        match self
        {
            TrussError::Number(number) => format!("Truss element with number {number} already exists!"),
            TrussError::SameNodes(node_1_number, node_2_number) => 
            {
                format!("Truss element with node number {node_1_number} and {node_2_number} already exists!")
            },
        }
    }
}


impl<V> FEM<V>
    where V: FloatTrait<Output = V>
{
    fn check_truss_data(&self, number: u32, node_1_number: u32, node_2_number: u32) -> Option<TrussError>
    {
        for (truss_number, truss_element) in self.get_truss_elements().iter()
        {
            if *truss_number == number
            {
                return Some(TrussError::Number(number));
            }
            if truss_element.is_nodes_numbers_same(node_1_number, node_2_number)
            {
                return Some(TrussError::SameNodes(node_1_number, node_2_number));
            }
        }
        None
    }

    
    pub fn add_truss(&mut self,
        number: u32,
        node_1_number: u32,
        node_2_number: u32,
        young_modulus: V,
        area: V,
        optional_area_2: Option<V>,
    )
        -> Result<(), String>
    {
        self.check_node_exist(node_1_number)?;
        self.check_node_exist(node_2_number)?;
            if let Some(truss_error) = self.check_truss_data(number, node_1_number, node_2_number)
        {
            return Err(truss_error.compose_error_message());
        }

        let truss_element = Truss::create(
            node_1_number, 
            node_2_number, 
            young_modulus, 
            area, 
            optional_area_2, 
            self.get_nodes(),
            self.get_props().get_rel_tol(),
            self.get_props().get_abs_tol(),
        )?;

        let rotation_matrix = truss_element.extract_rotation_matrix();
        let local_stiffness_matrix = truss_element.extract_local_stiffness_matrix(self.get_nodes())?;
        let transformed_local_stiffness_matrix = rotation_matrix
            .transpose()
            .multiply(&local_stiffness_matrix)?
            .multiply(&rotation_matrix)?;

        let node_1_index = self.get_nodes()
            .get(&node_1_number)
            .ok_or(format!("Node {node_1_number} is absent!"))?
            .get_index();
        let node_2_index = self.get_nodes()
            .get(&node_2_number)
            .ok_or(format!("Node {node_2_number} is absent!"))?
            .get_index();
        let start_positions = [
            ((0, 0), (node_1_index * NODE_DOF, node_1_index * NODE_DOF)),
            ((0, TRUSS_NODE_DOF), (node_1_index * NODE_DOF, node_2_index * NODE_DOF)),
            ((TRUSS_NODE_DOF, 0), (node_2_index * NODE_DOF, node_1_index * NODE_DOF)),
            ((TRUSS_NODE_DOF, TRUSS_NODE_DOF), (node_2_index * NODE_DOF, node_2_index * NODE_DOF)),
        ];

        for ((local_row, local_column), (global_row, global_column)) in start_positions
        {
            for i in 0..TRUSS_NODE_DOF
            {
                for j in 0..TRUSS_NODE_DOF
                {
                    let local_stiffness_matrix_element_value = transformed_local_stiffness_matrix
                        .get_element_value(&Position(local_row + i, local_column + j))?;
                    if *local_stiffness_matrix_element_value != V::from(0f32)
                    {
                        *self.get_mut_stiffness_matrix()
                            .get_mut_element_value(&Position(global_row + i, global_column + j))? += 
                                *local_stiffness_matrix_element_value;

                        if global_column + j >= global_row + i
                        {
                            let m = (global_column + j) - (global_row + i);
                            if m > self.get_skyline()[global_column + j]
                            {
                                self.get_mut_skyline()[global_column + j] = m;
                            }
                        }
                    }
                }
            }
        }

        self.get_mut_truss_elements().insert(number, truss_element);

        Ok(())
    }
}
