use extended_matrix::{BasicOperationsTrait, FloatTrait, Position};

use crate::fem::FEM;
use crate::fem::structs::{BEAM_NODE_DOF, Beam, NODE_DOF};

enum BeamElementError {
    Number(u32),
    SameNodes(u32, u32),
    NumberNotExist(u32),
}

impl BeamElementError {
    fn compose_error_message(&self) -> String {
        match self {
            Self::Number(number) => format!("Beam element with number {number} already exists!"),
            Self::SameNodes(node_1_number, node_2_number) => {
                format!(
                    "Beam element with node number {node_1_number} and {node_2_number} already exists!"
                )
            }
            Self::NumberNotExist(number) => {
                format!("Beam element with number {number} does not exist!")
            }
        }
    }
}

impl<V> FEM<V>
where
    V: FloatTrait<Output = V>,
{
    fn check_beam_data(
        &self,
        number: u32,
        node_1_number: u32,
        node_2_number: u32,
    ) -> Option<BeamElementError> {
        for (beam_number, beam_element) in self.get_beam_elements().iter() {
            if *beam_number == number {
                return Some(BeamElementError::Number(number));
            }
            if beam_element.is_nodes_numbers_same(node_1_number, node_2_number) {
                return Some(BeamElementError::SameNodes(node_1_number, node_2_number));
            }
        }
        None
    }

    pub fn add_beam(
        &mut self,
        number: u32,
        node_1_number: u32,
        node_2_number: u32,
        young_modulus: V,
        poisson_ratio: V,
        area: V,
        i11: V,
        i22: V,
        i12: V,
        it: V,
        shear_factor: V,
        local_axis_1_direction: [V; 3],
    ) -> Result<(), String> {
        self.check_node_exist(node_1_number)?;
        self.check_node_exist(node_2_number)?;
        if let Some(beam_error) = self.check_beam_data(number, node_1_number, node_2_number) {
            return Err(beam_error.compose_error_message());
        }

        let beam_element = Beam::create(
            number,
            node_1_number,
            node_2_number,
            young_modulus,
            poisson_ratio,
            area,
            i11,
            i22,
            i12,
            it,
            shear_factor,
            local_axis_1_direction,
            self.get_nodes(),
            self.get_props().get_rel_tol(),
            self.get_props().get_abs_tol(),
        )?;

        let rotation_matrix = beam_element.extract_rotation_matrix();

        let local_stiffness_matrix =
            beam_element.extract_local_stiffness_matrix(self.get_nodes())?;
        let transformed_local_stiffness_matrix = rotation_matrix
            .transpose()
            .multiply(&local_stiffness_matrix)?
            .multiply(&rotation_matrix)?;

        let node_1_index = self
            .get_nodes()
            .get(&node_1_number)
            .ok_or(format!("Node {node_1_number} is absent!"))?
            .get_index();
        let node_2_index = self
            .get_nodes()
            .get(&node_2_number)
            .ok_or(format!("Node {node_2_number} is absent!"))?
            .get_index();
        let start_positions = [
            ((0, 0), (node_1_index * NODE_DOF, node_1_index * NODE_DOF)),
            (
                (0, BEAM_NODE_DOF),
                (node_1_index * NODE_DOF, node_2_index * NODE_DOF),
            ),
            (
                (BEAM_NODE_DOF, 0),
                (node_2_index * NODE_DOF, node_1_index * NODE_DOF),
            ),
            (
                (BEAM_NODE_DOF, BEAM_NODE_DOF),
                (node_2_index * NODE_DOF, node_2_index * NODE_DOF),
            ),
        ];

        for ((local_row, local_column), (global_row, global_column)) in start_positions {
            for i in 0..BEAM_NODE_DOF {
                for j in 0..BEAM_NODE_DOF {
                    let local_stiffness_matrix_element_value =
                        transformed_local_stiffness_matrix
                            .get_element_value(&Position(local_row + i, local_column + j))?;
                    if *local_stiffness_matrix_element_value != V::from(0f32) {
                        *self
                            .get_mut_stiffness_matrix()
                            .get_mut_element_value(&Position(
                                global_row + i,
                                global_column + j,
                            ))? += *local_stiffness_matrix_element_value;
                    }
                }
            }
        }

        self.get_mut_beam_elements().insert(number, beam_element);

        Ok(())
    }

    pub(crate) fn check_beam_element_exist(&self, number: u32) -> Result<(), String> {
        if !self.get_beam_elements().contains_key(&number) {
            return Err(BeamElementError::NumberNotExist(number).compose_error_message());
        }
        Ok(())
    }
}
