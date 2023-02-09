use extended_matrix::{FloatTrait, BasicOperationsTrait, Position, Matrix, SquareMatrix};

use crate::fem::FEM;
use crate::fem::structs::{NODE_DOF, SeparatedStiffnessMatrix};
use crate::fem::methods_for_bc_data_handle::DOFParameter;


enum SeparatingMatrixError
{
    Load(u32, DOFParameter),
    Displacement(u32, DOFParameter),
    NoRestraints,
}


impl SeparatingMatrixError
{
    fn compose_error_message(&self) -> String
    {
        match self
        {
            Self::Load(node_number, dof_parameter) => 
            {
                format!("There are no stiffness to withstand load {dof_parameter:?} applied to node {node_number}!")
            },
            Self::Displacement(node_number, dof_parameter) => 
            {
                format!("There are no stiffness to withstand displacement {dof_parameter:?} applied \
                    to node {node_number}!")
            },
            Self::NoRestraints => format!("There are no restraints applied!")
        }
    }
}


impl<V> FEM<V>
    where V: FloatTrait<Output = V>
{
    fn check_excluded_index(&self, index: usize) -> Result<(), String>
    {
        let node_number = self.get_index_node_number_map()
            .get(&(index / NODE_DOF))
            .ok_or(format!("Node with index {index} is absent!"))?;
        let dof_parameter = DOFParameter::from_usize(index % NODE_DOF);

        if self.get_imposed_constraints()[index]
        {
            return Err(SeparatingMatrixError::Displacement(*node_number, dof_parameter).compose_error_message());
        }

        if *self.get_forces_vector().get_element_value(&Position(index, 0))? != V::from(0f32)
        {
            return Err(SeparatingMatrixError::Load(*node_number, dof_parameter).compose_error_message());
        }

        Ok(())
    }


    pub(crate) fn separate_stiffness_matrix(&self) -> Result<SeparatedStiffnessMatrix<V>, String>
    {
        let mut k_aa_indexes = Vec::new();
        let mut k_bb_indexes = Vec::new();

        for index in self.get_indexes().iter()
        {
            if *self.get_stiffness_matrix().get_element_value(&Position(*index, *index))? == V::from(0f32)
            {
                self.check_excluded_index(*index)?;
            }
            else
            {
                if self.get_imposed_constraints()[*index]
                {
                    k_bb_indexes.push(*index);
                }
                else
                {
                    k_aa_indexes.push(*index);
                }
            }
        }

        if k_bb_indexes.len() == 0
        {
            return Err(SeparatingMatrixError::NoRestraints.compose_error_message());
        }

        let mut k_aa_skyline = vec![0; k_aa_indexes.len()];

        let mut k_aa_matrix = SquareMatrix::create(
            k_aa_indexes.len(), 
            &vec![V::from(0f32); k_aa_indexes.len() * k_aa_indexes.len()],
        );
        let mut k_ab_matrix = Matrix::create(
            k_aa_indexes.len(),
            k_bb_indexes.len(),
            &vec![V::from(0f32); k_aa_indexes.len() * k_bb_indexes.len()],
        );
        let mut k_ba_matrix = Matrix::create(
            k_bb_indexes.len(),
            k_aa_indexes.len(),
            &vec![V::from(0f32); k_bb_indexes.len() * k_aa_indexes.len()],
        );
        let mut k_bb_matrix = SquareMatrix::create(
            k_bb_indexes.len(), 
            &vec![V::from(0f32); k_bb_indexes.len() * k_bb_indexes.len()],
        );
        for i in 0..k_aa_indexes.len()
        {
            let matrix_row = k_aa_indexes[i];
            let matrix_column = k_aa_indexes[i];

            *k_aa_matrix.get_mut_element_value(&Position(i, i))? = *self.get_stiffness_matrix()
                .get_element_value(&Position(matrix_row, matrix_column))?;

            for j in i + 1..k_aa_indexes.len()
            {
                let matrix_row = k_aa_indexes[i];
                let matrix_column = k_aa_indexes[j];
                *k_aa_matrix.get_mut_element_value(&Position(i, j))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(matrix_row, matrix_column))?;
                *k_aa_matrix.get_mut_element_value(&Position(j, i))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(matrix_column, matrix_row))?;

                if *self.get_stiffness_matrix()
                    .get_element_value(&Position(matrix_row, matrix_column))? != V::from(0f32)
                {
                    if j >= i && j - i > k_aa_skyline[j]
                    {
                        k_aa_skyline[j] = j - i;
                    }
                }
            }
        }

        for i in 0..k_aa_indexes.len()
        {
            let mut row = k_aa_indexes[i];
            let mut column = k_aa_indexes[i];

            *k_aa_matrix.get_mut_element_value(&Position(i, i))? = *self.get_stiffness_matrix()
                .get_element_value(&Position(row, column))?;

            for j in i + 1..k_aa_indexes.len()
            {
                row = k_aa_indexes[i];
                column = k_aa_indexes[j];
                *k_aa_matrix.get_mut_element_value(&Position(i, j))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(row, column))?;
                *k_aa_matrix.get_mut_element_value(&Position(j, i))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(column, row))?;

                if *self.get_stiffness_matrix()
                    .get_element_value(&Position(row, column))? != V::from(0f32)
                {
                    if j >= i && j - i > k_aa_skyline[j]
                    {
                        k_aa_skyline[j] = j - i;
                    }
                }
            }

            for k in 0..k_bb_indexes.len()
            {
                row = k_aa_indexes[i];
                column = k_bb_indexes[k];

                *k_ab_matrix.get_mut_element_value(&Position(i, k))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(row, column))?;
                *k_ba_matrix.get_mut_element_value(&Position(k, i))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(column, row))?;

                row = k_bb_indexes[k];
                column = k_bb_indexes[k];

                *k_bb_matrix.get_mut_element_value(&Position(k, k))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(row, column))?;

                for l in k + 1..k_bb_indexes.len()
                {
                    row = k_bb_indexes[k];
                    column = k_bb_indexes[l];

                    *k_bb_matrix.get_mut_element_value(&Position(k, l))? = *self.get_stiffness_matrix()
                        .get_element_value(&Position(row, column))?;
                    *k_bb_matrix.get_mut_element_value(&Position(l, k))? = *self.get_stiffness_matrix()
                        .get_element_value(&Position(column, row))?;
                }
            }
        }
        Ok(
            SeparatedStiffnessMatrix::create(
                k_aa_indexes, k_bb_indexes, k_aa_skyline, k_aa_matrix, k_ab_matrix, k_ba_matrix, k_bb_matrix,
            )
        )
    }
}
