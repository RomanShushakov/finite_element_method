use extended_matrix::{FloatTrait, BasicOperationsTrait, Position};

use crate::fem::FEM;
use crate::fem::structs::{Plate, PLATE_NODE_DOF, NODE_DOF};


enum PlateElementError
{
    Number(u32),
    SameNodes(u32, u32, u32, u32),
    NumberNotExist(u32),
}


impl PlateElementError
{
    fn compose_error_message(&self) -> String
    {
        match self
        {
            Self::Number(number) => format!("Plate element with number {number} already exists!"),
            Self::SameNodes(node_1_number, node_2_number, node_3_number, node_4_number) => 
            {
                format!("Plate element with nodes numbers [{node_1_number}, {node_2_number}, \
                    {node_3_number}, {node_4_number}] already exists!")
            },
            Self::NumberNotExist(number) => format!("Plate element with number {number} does not exist!"),
        }
    }
}


impl<V> FEM<V>
    where V: FloatTrait<Output = V>
{
    fn check_plate_data(
        &self, number: u32, node_1_number: u32, node_2_number: u32, node_3_number: u32, node_4_number: u32,
    ) 
        -> Option<PlateElementError>
    {
        for (plate_number, plate_element) in self.get_plate_elements().iter()
        {
            if *plate_number == number
            {
                return Some(PlateElementError::Number(number));
            }
            if plate_element.is_nodes_numbers_same(
                &[node_1_number, node_2_number, node_3_number, node_4_number],
            )
            {
                return Some(PlateElementError::SameNodes(node_1_number, node_2_number, node_3_number, node_4_number));
            }
        }
        None
    }


    pub fn add_plate(
        &mut self,
        number: u32,
        node_1_number: u32,
        node_2_number: u32,
        node_3_number: u32,
        node_4_number: u32,
        young_modulus: V,
        poisson_ratio: V,
        thickness: V,
        shear_factor: V,
    )
        -> Result<(), String>
    {
        self.check_node_exist(node_1_number)?;
        self.check_node_exist(node_2_number)?;
        self.check_node_exist(node_3_number)?;
        self.check_node_exist(node_4_number)?;

        if let Some(plate_error) = self.check_plate_data(
            number, node_1_number, node_2_number, node_3_number, node_4_number,
        )
        {
            return Err(plate_error.compose_error_message());
        }

        let plate_element = Plate::create(
            number,
            node_1_number,
            node_2_number,
            node_3_number,
            node_4_number,
            young_modulus,
            poisson_ratio,
            thickness,
            shear_factor,
            self.get_nodes(),
            self.get_props().get_rel_tol(),
            self.get_props().get_abs_tol(),
        )?;

        let rotation_matrix = plate_element.extract_rotation_matrix();

        let f = |data: &str| println!("{data}");
        rotation_matrix.show(f);

        // let local_stiffness_matrix = beam_element.extract_local_stiffness_matrix(self.get_nodes())?;
        // let transformed_local_stiffness_matrix = rotation_matrix
        //     .transpose()
        //     .multiply(&local_stiffness_matrix)?
        //     .multiply(&rotation_matrix)?;

        // let node_1_index = self.get_nodes()
        //     .get(&node_1_number)
        //     .ok_or(format!("Node {node_1_number} is absent!"))?
        //     .get_index();
        // let node_2_index = self.get_nodes()
        //     .get(&node_2_number)
        //     .ok_or(format!("Node {node_2_number} is absent!"))?
        //     .get_index();
        // let start_positions = [
        //     ((0, 0), (node_1_index * NODE_DOF, node_1_index * NODE_DOF)),
        //     ((0, BEAM_NODE_DOF), (node_1_index * NODE_DOF, node_2_index * NODE_DOF)),
        //     ((BEAM_NODE_DOF, 0), (node_2_index * NODE_DOF, node_1_index * NODE_DOF)),
        //     ((BEAM_NODE_DOF, BEAM_NODE_DOF), (node_2_index * NODE_DOF, node_2_index * NODE_DOF)),
        // ];

        // for ((local_row, local_column), (global_row, global_column)) in start_positions
        // {
        //     for i in 0..BEAM_NODE_DOF
        //     {
        //         for j in 0..BEAM_NODE_DOF
        //         {
        //             let local_stiffness_matrix_element_value = transformed_local_stiffness_matrix
        //                 .get_element_value(&Position(local_row + i, local_column + j))?;
        //             if *local_stiffness_matrix_element_value != V::from(0f32)
        //             {
        //                 *self.get_mut_stiffness_matrix()
        //                     .get_mut_element_value(&Position(global_row + i, global_column + j))? += 
        //                         *local_stiffness_matrix_element_value;
        //             }
        //         }
        //     }
        // }

        // self.get_mut_beam_elements().insert(number, beam_element);

        Ok(())
    }


    // pub(crate) fn check_beam_element_exist(&self, number: u32) -> Result<(), String>
    // {
    //     if !self.get_beam_elements().contains_key(&number) 
    //     {
    //         return Err(BeamElementError::NumberNotExist(number).compose_error_message());
    //     }
    //     Ok(())
    // }
}
