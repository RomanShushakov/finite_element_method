use extended_matrix::{BasicOperationsTrait, FloatTrait, Matrix, Position, SquareMatrix, Vector};

use crate::fem::FEM;
use crate::fem::methods_for_bc_data_handle::DOFParameter;
use crate::fem::structs::{NODE_DOF, SeparatedStiffnessMatrix, SeparatedStiffnessMatrixSparse};

enum SeparatingMatrixError {
    Load(u32, DOFParameter),
    Displacement(u32, DOFParameter),
    NoRestraints,
}

impl SeparatingMatrixError {
    fn compose_error_message(&self) -> String {
        match self {
            Self::Load(node_number, dof_parameter) => {
                format!(
                    "There are no stiffness to withstand load {dof_parameter:?} applied to node {node_number}!"
                )
            }
            Self::Displacement(node_number, dof_parameter) => {
                format!(
                    "There are no stiffness to withstand displacement {dof_parameter:?} applied \
                    to node {node_number}!"
                )
            }
            Self::NoRestraints => format!("There are no restraints applied!"),
        }
    }
}

impl<V> FEM<V>
where
    V: FloatTrait<Output = V>,
{
    fn check_excluded_index(&self, index: usize) -> Result<(), String> {
        let node_number = self
            .get_index_node_number_map()
            .get(&(index / NODE_DOF))
            .ok_or(format!("Node with index {index} is absent!"))?;
        let dof_parameter = DOFParameter::from_usize(index % NODE_DOF);

        if self.get_imposed_constraints()[index] {
            return Err(
                SeparatingMatrixError::Displacement(*node_number, dof_parameter)
                    .compose_error_message(),
            );
        }

        if *self
            .get_forces_vector()
            .get_element_value(&Position(index, 0))?
            != V::from(0f32)
        {
            return Err(
                SeparatingMatrixError::Load(*node_number, dof_parameter).compose_error_message()
            );
        }

        Ok(())
    }

    pub fn separate_stiffness_matrix_direct(
        &mut self,
    ) -> Result<SeparatedStiffnessMatrix<V>, String> {
        let dense_stiffness_matrix = self.get_stiffness_matrix().to_dense();
        *self.get_mut_stiffness_matrix() = dense_stiffness_matrix;

        let mut k_aa_indexes = Vec::new();
        let mut k_bb_indexes = Vec::new();

        for index in self.get_indexes().iter() {
            if *self
                .get_stiffness_matrix()
                .get_element_value(&Position(*index, *index))?
                == V::from(0f32)
            {
                self.check_excluded_index(*index)?;
            } else {
                if self.get_imposed_constraints()[*index] {
                    k_bb_indexes.push(*index);
                } else {
                    k_aa_indexes.push(*index);
                }
            }
        }

        if k_bb_indexes.len() == 0 {
            return Err(SeparatingMatrixError::NoRestraints.compose_error_message());
        }

        let mut k_aa_skyline = vec![0; k_aa_indexes.len()];

        let mut k_aa_matrix = SquareMatrix::create(1, &[V::from(0f32)]);
        k_aa_matrix
            .get_mut_shape()
            .update(k_aa_indexes.len(), k_aa_indexes.len());

        let mut k_ab_matrix = Matrix::create(1, 1, &[V::from(0f32)]);
        k_ab_matrix
            .get_mut_shape()
            .update(k_aa_indexes.len(), k_bb_indexes.len());

        let mut k_ba_matrix = Matrix::create(1, 1, &[V::from(0f32)]);
        k_ba_matrix
            .get_mut_shape()
            .update(k_bb_indexes.len(), k_aa_indexes.len());

        let mut k_bb_matrix = SquareMatrix::create(1, &[V::from(0f32)]);
        k_bb_matrix
            .get_mut_shape()
            .update(k_bb_indexes.len(), k_bb_indexes.len());

        for i in 0..k_aa_indexes.len() {
            let mut row = k_aa_indexes[i];
            let mut column = k_aa_indexes[i];

            let mut value = self
                .get_mut_stiffness_matrix()
                .get_mut_elements()
                .remove(&Position(row, column))
                .ok_or(format!("Element ({row}, {column}) is absent!"))?;
            k_aa_matrix.get_mut_elements().insert(Position(i, i), value);

            for j in i + 1..k_aa_indexes.len() {
                row = k_aa_indexes[i];
                column = k_aa_indexes[j];

                if *self
                    .get_stiffness_matrix()
                    .get_element_value(&Position(row, column))?
                    != V::from(0f32)
                {
                    if j >= i && j - i > k_aa_skyline[j] {
                        k_aa_skyline[j] = j - i;
                    }
                }
                value = self
                    .get_mut_stiffness_matrix()
                    .get_mut_elements()
                    .remove(&Position(row, column))
                    .ok_or(format!("Element ({row}, {column}) is absent!"))?;
                k_aa_matrix.get_mut_elements().insert(Position(i, j), value);

                value = self
                    .get_mut_stiffness_matrix()
                    .get_mut_elements()
                    .remove(&Position(column, row))
                    .ok_or(format!("Element ({column}, {row}) is absent!"))?;
                k_aa_matrix.get_mut_elements().insert(Position(j, i), value);
            }

            for k in 0..k_bb_indexes.len() {
                row = k_aa_indexes[i];
                column = k_bb_indexes[k];

                value = self
                    .get_mut_stiffness_matrix()
                    .get_mut_elements()
                    .remove(&Position(row, column))
                    .ok_or(format!("Element ({row}, {column}) is absent!"))?;
                k_ab_matrix.get_mut_elements().insert(Position(i, k), value);

                value = self
                    .get_mut_stiffness_matrix()
                    .get_mut_elements()
                    .remove(&Position(column, row))
                    .ok_or(format!("Element ({column}, {row}) is absent!"))?;
                k_ba_matrix.get_mut_elements().insert(Position(k, i), value);
            }
        }

        for k in 0..k_bb_indexes.len() {
            let mut row = k_bb_indexes[k];
            let mut column = k_bb_indexes[k];

            let mut value = self
                .get_mut_stiffness_matrix()
                .get_mut_elements()
                .remove(&Position(row, column))
                .ok_or(format!("Element ({row}, {column}) is absent!"))?;
            k_bb_matrix.get_mut_elements().insert(Position(k, k), value);

            for l in k + 1..k_bb_indexes.len() {
                row = k_bb_indexes[k];
                column = k_bb_indexes[l];

                value = self
                    .get_mut_stiffness_matrix()
                    .get_mut_elements()
                    .remove(&Position(row, column))
                    .ok_or(format!("Element ({row}, {column}) is absent!"))?;
                k_bb_matrix.get_mut_elements().insert(Position(k, l), value);

                value = self
                    .get_mut_stiffness_matrix()
                    .get_mut_elements()
                    .remove(&Position(column, row))
                    .ok_or(format!("Element ({column}, {row}) is absent!"))?;
                k_bb_matrix.get_mut_elements().insert(Position(l, k), value);
            }
        }

        *self.get_mut_stiffness_matrix() = SquareMatrix::create(1, &[V::from(0f32)]);

        Ok(SeparatedStiffnessMatrix::create(
            k_aa_indexes,
            k_bb_indexes,
            k_aa_skyline,
            k_aa_matrix,
            k_ab_matrix,
            k_ba_matrix,
            k_bb_matrix,
        ))
    }

    pub fn separate_stiffness_matrix_sparse_iterative(
        &self,
    ) -> Result<SeparatedStiffnessMatrixSparse<V>, String> {
        let n = self.get_indexes().len();

        let mut diag = vec![V::from(0.0_f32); n];
        for (pos, val) in self.get_stiffness_matrix().get_elements().iter() {
            if pos.0 == pos.1 {
                diag[pos.0] = *val;
            }
        }

        let mut k_aa_indexes = Vec::new();
        let mut k_bb_indexes = Vec::new();

        for &idx in self.get_indexes().iter() {
            if diag[idx] == V::from(0.0_f32) {
                if self.get_imposed_constraints()[idx] {
                    return Err(format!(
                        "There are no stiffness to withstand displacement {:?} applied to node {}!",
                        DOFParameter::from_usize(idx % NODE_DOF),
                        self.get_index_node_number_map()
                            .get(&(idx / NODE_DOF))
                            .copied()
                            .unwrap_or(0)
                    ));
                } else {
                    // excluded/inactive dof
                    continue;
                }
            }

            if self.get_imposed_constraints()[idx] {
                k_bb_indexes.push(idx);
            } else {
                k_aa_indexes.push(idx);
            }
        }

        if k_bb_indexes.is_empty() {
            return Err("No restraints".to_string());
        }

        let n_aa = k_aa_indexes.len();
        let n_bb = k_bb_indexes.len();

        let mut aa_pos = vec![usize::MAX; n];
        let mut bb_pos = vec![usize::MAX; n];
        for (li, &gi) in k_aa_indexes.iter().enumerate() {
            aa_pos[gi] = li;
        }
        for (lj, &gj) in k_bb_indexes.iter().enumerate() {
            bb_pos[gj] = lj;
        }

        let mut k_aa_triplets = Vec::new();
        let mut k_ab_triplets = Vec::new();
        let mut k_ba_triplets = Vec::new();
        let mut k_bb_triplets = Vec::new();

        for (pos, val) in self.get_stiffness_matrix().get_elements().iter() {
            if *val == V::from(0.0_f32) {
                continue;
            }

            let gi = pos.0;
            let gj = pos.1;

            let i_aa = aa_pos[gi];
            let j_aa = aa_pos[gj];
            let i_bb = bb_pos[gi];
            let j_bb = bb_pos[gj];

            if i_aa != usize::MAX && j_aa != usize::MAX {
                k_aa_triplets.push((i_aa, j_aa, *val));
            } else if i_aa != usize::MAX && j_bb != usize::MAX {
                k_ab_triplets.push((i_aa, j_bb, *val));
            } else if i_bb != usize::MAX && j_aa != usize::MAX {
                k_ba_triplets.push((i_bb, j_aa, *val));
            } else if i_bb != usize::MAX && j_bb != usize::MAX {
                k_bb_triplets.push((i_bb, j_bb, *val));
            } else {
                // entry touches excluded dof(s) -> ignore
                continue;
            }
        }

        if k_aa_triplets.is_empty() {
            return Err(
                "Sparse separation: K_aa is empty (structure has no free stiffness?)".to_string(),
            );
        }

        Ok(SeparatedStiffnessMatrixSparse::create(
            k_aa_indexes,
            k_bb_indexes,
            k_aa_triplets,
            k_ab_triplets,
            k_ba_triplets,
            k_bb_triplets,
            n_aa,
            n_bb,
        ))
    }

    pub fn compose_r_a_vector(&self, k_aa_indexes: &Vec<usize>) -> Result<Vector<V>, String> {
        let mut r_a_vector = Vector::create(&vec![V::from(0f32); k_aa_indexes.len()]);
        for i in 0..k_aa_indexes.len() {
            *r_a_vector.get_mut_element_value(&Position(i, 0))? = *self
                .get_forces_vector()
                .get_element_value(&Position(k_aa_indexes[i], 0))?;
        }

        Ok(r_a_vector)
    }

    pub fn compose_u_b_vector(&self, k_bb_indexes: &Vec<usize>) -> Result<Vector<V>, String> {
        let mut u_b_vector = Vector::create(&vec![V::from(0f32); k_bb_indexes.len()]);
        for i in 0..k_bb_indexes.len() {
            *u_b_vector.get_mut_element_value(&Position(i, 0))? = *self
                .get_displacements_vector()
                .get_element_value(&Position(k_bb_indexes[i], 0))?;
        }

        Ok(u_b_vector)
    }
}
