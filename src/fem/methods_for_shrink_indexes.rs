use extended_matrix::{FloatTrait, BasicOperationsTrait, Position};

use crate::fem::FEM;
// use crate::fem::structs::NODE_DOF;


// #[derive(Debug, Clone, Copy)]
// pub enum DOFParameter
// {
//     X = 0,
//     Y = 1,
//     Z = 2,
//     ThX = 3,
//     ThY = 4,
//     ThZ = 5,
// }


impl<V> FEM<V>
    where V: FloatTrait<Output = V>
{
    pub(crate) fn shrink_indexes(&self) -> Result<(), String>
    {
        
        Ok(())
    }
    // pub fn add_concentrated_load(
    //     &mut self, node_number: u32, dof_parameter: DOFParameter, value: V,
    // )
    //     -> Result<(), String>
    // {
    //     self.check_node_exist(node_number)?;

    //     let node_index = self.get_nodes()
    //         .get(&node_number)
    //         .ok_or(format!("Node {node_number} is absent!"))?
    //         .get_index();

    //     *self.get_mut_forces_vector()
    //         .get_mut_element_value(&Position(node_index * NODE_DOF + dof_parameter as usize, 0))? += value;

    //     Ok(())
    // }


    // pub fn add_displacement(
    //     &mut self, node_number: u32, dof_parameter: DOFParameter, value: V,
    // )
    //     -> Result<(), String>
    // {
    //     self.check_node_exist(node_number)?;

    //     let node_index = self.get_nodes()
    //         .get(&node_number)
    //         .ok_or(format!("Node {node_number} is absent!"))?
    //         .get_index();

    //     let index = node_index * NODE_DOF + dof_parameter as usize;
    //     if self.get_imposed_constraints()[index]
    //     {
    //         return Err(format!("Displacement {dof_parameter:?} already applied to node {node_number}!"));
    //     }

    //     self.get_mut_imposed_constraints()[index] = true;
    //     *self.get_mut_displacements_vector()
    //         .get_mut_element_value(&Position(index, 0))? = value;

    //     Ok(())
    // }
}
