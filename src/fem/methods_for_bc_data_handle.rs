use extended_matrix::{FloatTrait, BasicOperationsTrait, Position};

use crate::fem::FEM;
use crate::fem::structs::NODE_DOF;


#[derive(Debug, Clone, Copy)]
pub enum DOFParameter
{
    X = 0,
    Y = 1,
    Z = 2,
    ThX = 3,
    ThY = 4,
    ThZ = 5,
}


impl DOFParameter
{
    pub fn from_usize(number: usize) -> Self
    {
        match number
        {
            0 => DOFParameter::X,
            1 => DOFParameter::Y,
            2 => DOFParameter::Z,
            3 => DOFParameter::ThX,
            4 => DOFParameter::ThY,
            5 => DOFParameter::ThZ,
            _ => unreachable!("Incorrect dof parameter number!"),
        }
    }
}


impl<V> FEM<V>
    where V: FloatTrait<Output = V>
{
    pub fn add_concentrated_load(
        &mut self, node_number: u32, dof_parameter: DOFParameter, value: V,
    )
        -> Result<(), String>
    {
        self.check_node_exist(node_number)?;

        let node_index = self.get_nodes()
            .get(&node_number)
            .ok_or(format!("Node {node_number} is absent!"))?
            .get_index();

        *self.get_mut_forces_vector()
            .get_mut_element_value(&Position(node_index * NODE_DOF + dof_parameter as usize, 0))? += value;

        Ok(())
    }


    pub fn add_uniformly_distributed_line_load(
        &mut self, beam_element_number: u32, dof_parameter: DOFParameter, value: V,
    )
        -> Result<(), String>
    {
        self.check_beam_element_exist(beam_element_number)?;

        let beam = self.get_beam_elements()
            .get(&beam_element_number)
            .ok_or(format!("Beam element {beam_element_number} is absent!"))?;

        let nodal_forces = beam.convert_uniformly_distributed_line_load_to_nodal_loads(
            value, &self.get_nodes(),
        )?;

        let [node_1_number, node_2_number] = beam.get_nodes_numbers();

        let node_1_index = self.get_nodes()
            .get(&node_1_number)
            .ok_or(format!("Node {node_1_number} is absent!"))?
            .get_index();
        *self.get_mut_forces_vector()
            .get_mut_element_value(&Position(node_1_index * NODE_DOF + dof_parameter as usize, 0))? += 
                *nodal_forces.get_element_value(&Position(0, 0))?;

        let node_2_index = self.get_nodes()
            .get(&node_2_number)
            .ok_or(format!("Node {node_2_number} is absent!"))?
            .get_index();
        *self.get_mut_forces_vector()
            .get_mut_element_value(&Position(node_2_index * NODE_DOF + dof_parameter as usize, 0))? += 
                *nodal_forces.get_element_value(&Position(1, 0))?;

        Ok(())
    }


    pub fn add_uniformly_distributed_surface_load(
        &mut self, plate_element_number: u32, dof_parameter: DOFParameter, value: V,
    )
        -> Result<(), String>
    {
        self.check_plate_element_exist(plate_element_number)?;

        let plate = self.get_plate_elements()
            .get(&plate_element_number)
            .ok_or(format!("Plate element {plate_element_number} is absent!"))?;

        let nodal_forces = plate.convert_uniformly_distributed_surface_load_to_nodal_loads(
            value, &self.get_nodes(), self.get_props().get_rel_tol(),
        )?;

        let [node_1_number, node_2_number, node_3_number, node_4_number] = plate.get_nodes_numbers();

        let node_1_index = self.get_nodes()
            .get(&node_1_number)
            .ok_or(format!("Node {node_1_number} is absent!"))?
            .get_index();
        *self.get_mut_forces_vector()
            .get_mut_element_value(&Position(node_1_index * NODE_DOF + dof_parameter as usize, 0))? += 
                *nodal_forces.get_element_value(&Position(0, 0))?;

        let node_2_index = self.get_nodes()
            .get(&node_2_number)
            .ok_or(format!("Node {node_2_number} is absent!"))?
            .get_index();
        *self.get_mut_forces_vector()
            .get_mut_element_value(&Position(node_2_index * NODE_DOF + dof_parameter as usize, 0))? += 
                *nodal_forces.get_element_value(&Position(1, 0))?;
            
        let node_3_index = self.get_nodes()
            .get(&node_3_number)
            .ok_or(format!("Node {node_3_number} is absent!"))?
            .get_index();
        *self.get_mut_forces_vector()
            .get_mut_element_value(&Position(node_3_index * NODE_DOF + dof_parameter as usize, 0))? += 
                *nodal_forces.get_element_value(&Position(2, 0))?;

        let node_4_index = self.get_nodes()
            .get(&node_4_number)
            .ok_or(format!("Node {node_4_number} is absent!"))?
            .get_index();
        *self.get_mut_forces_vector()
            .get_mut_element_value(&Position(node_4_index * NODE_DOF + dof_parameter as usize, 0))? += 
                *nodal_forces.get_element_value(&Position(3, 0))?;

        Ok(())
    }


    pub fn add_displacement(
        &mut self, node_number: u32, dof_parameter: DOFParameter, value: V,
    )
        -> Result<(), String>
    {
        self.check_node_exist(node_number)?;

        let node_index = self.get_nodes()
            .get(&node_number)
            .ok_or(format!("Node {node_number} is absent!"))?
            .get_index();

        let index = node_index * NODE_DOF + dof_parameter as usize;
        if self.get_imposed_constraints()[index]
        {
            return Err(format!("Displacement {dof_parameter:?} already applied to node {node_number}!"));
        }

        self.get_mut_imposed_constraints()[index] = true;
        *self.get_mut_displacements_vector()
            .get_mut_element_value(&Position(index, 0))? = value;

        Ok(())
    }
}
