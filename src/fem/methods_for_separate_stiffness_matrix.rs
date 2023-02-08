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
        // let mut shrinked_indexes = Vec::new();
        let mut k_aa_indexes = Vec::new();
        let mut k_bb_indexes = Vec::new();
        // let mut shrinked_imposed_constraints = Vec::new();
        let mut skyline = self.get_skyline().clone();
        // let mut interim_kaa_skyline = self.get_skyline().clone();
        // let mut shrinked_skyline = Vec::new();
        let mut k_aa_skyline = Vec::new();
        for index in self.get_indexes().iter()
        {
            if *self.get_stiffness_matrix().get_element_value(&Position(*index, *index))? == V::from(0f32)
            {
                self.check_excluded_index(*index)?;
                for i in *index..skyline.len()
                {
                    if skyline[i] != 0
                    {
                        if *index >= i - skyline[i]
                        {
                            skyline[i] -= 1;
                        }
                    }
                }
            }
            else
            {
                // shrinked_indexes.push(*index);
                // shrinked_imposed_constraints.push(self.get_imposed_constraints()[*index]);
                // shrinked_skyline.push(skyline[*index]);

                if self.get_imposed_constraints()[*index]
                {
                    k_bb_indexes.push(*index);
                    for i in *index..skyline.len()
                    {
                        if skyline[i] != 0
                        {
                            if *index >= i - skyline[i]
                            {
                                skyline[i] -= 1;
                            }
                        }
                    }
                }
                else
                {
                    k_aa_indexes.push(*index);
                    k_aa_skyline.push(skyline[*index]);
                }
            }
        }


        // let mut shrinked_matrix = SquareMatrix::create(
        //     shrinked_indexes.len(), 
        //     &vec![V::from(0f32); shrinked_indexes.len() * shrinked_indexes.len()],
        // );
        // for i in 0..shrinked_indexes.len()
        // {
        //     let matrix_row = shrinked_indexes[i];
        //     let matrix_column = shrinked_indexes[i];

        //     *shrinked_matrix.get_mut_element_value(&Position(i, i))? = *self.get_stiffness_matrix()
        //         .get_element_value(&Position(matrix_row, matrix_column))?;

        //     for j in i + 1..shrinked_indexes.len()
        //     {
        //         let matrix_row = shrinked_indexes[i];
        //         let matrix_column = shrinked_indexes[j];
        //         *shrinked_matrix.get_mut_element_value(&Position(i, j))? = *self.get_stiffness_matrix()
        //             .get_element_value(&Position(matrix_row, matrix_column))?;
        //         *shrinked_matrix.get_mut_element_value(&Position(j, i))? = *self.get_stiffness_matrix()
        //             .get_element_value(&Position(matrix_column, matrix_row))?;
        //     }
        // }

        // println!("Shrinked indexes:\n{shrinked_indexes:?}\n");
        // println!("Shrinked imposed constraints:\n{shrinked_imposed_constraints:?}\n");
        // println!("Shrinked skyline:\n{shrinked_skyline:?}\n");
        // println!("Shrinked matrix:");
        // let f = |data: &str| println!("{data}");
        // shrinked_matrix.show(f);
        // println!();

        if k_bb_indexes.len() == 0
        {
            return Err(SeparatingMatrixError::NoRestraints.compose_error_message());
        }

        let mut k_aa_matrix = SquareMatrix::create(
            k_aa_indexes.len(), 
            &vec![V::from(0f32); k_aa_indexes.len() * k_aa_indexes.len()],
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
            }
        }

        for i in 0..k_aa_skyline.len()
        {
            let mut m = k_aa_skyline[i];
            while *k_aa_matrix.get_element_value(&Position(i - m, i))? == V::from(0f32) && m != 0
            {
                m -= 1;
            }
            k_aa_skyline[i] = m;
        }

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
        for i in 0..k_aa_indexes.len()
        {
            for j in 0..k_bb_indexes.len()
            {
                let matrix_row = k_aa_indexes[i];
                let matrix_column = k_bb_indexes[j];

                *k_ab_matrix.get_mut_element_value(&Position(i, j))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(matrix_row, matrix_column))?;
                *k_ba_matrix.get_mut_element_value(&Position(j, i))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(matrix_column, matrix_row))?;
            }
        }

        let mut k_bb_matrix = SquareMatrix::create(
            k_bb_indexes.len(), 
            &vec![V::from(0f32); k_bb_indexes.len() * k_bb_indexes.len()],
        );
        for i in 0..k_bb_indexes.len()
        {
            let matrix_row = k_bb_indexes[i];
            let matrix_column = k_bb_indexes[i];

            *k_bb_matrix.get_mut_element_value(&Position(i, i))? = *self.get_stiffness_matrix()
                .get_element_value(&Position(matrix_row, matrix_column))?;

            for j in i + 1..k_bb_indexes.len()
            {
                let matrix_row = k_bb_indexes[i];
                let matrix_column = k_bb_indexes[j];

                *k_bb_matrix.get_mut_element_value(&Position(i, j))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(matrix_row, matrix_column))?;
                *k_bb_matrix.get_mut_element_value(&Position(j, i))? = *self.get_stiffness_matrix()
                    .get_element_value(&Position(matrix_column, matrix_row))?;
            }
        }


        // println!("Kaa indexes:\n{k_aa_indexes:?}\n");
        // println!("Kbb indexes:\n{k_bb_indexes:?}\n");
        // println!("Kaa skyline:\n{k_aa_skyline:?}\n");
        // let f = |data: &str| println!("{data}");
        // println!("Kaa matrix:");
        // k_aa_matrix.show(f);
        // println!();
        // println!("Kab matrix:");
        // k_ab_matrix.show(f);
        // println!();
        // println!("Kba matrix:");
        // k_ba_matrix.show(f);
        // println!();
        // println!("Kbb matrix:");
        // k_bb_matrix.show(f);
        // println!();

        Ok(
            SeparatedStiffnessMatrix::create(
                k_aa_indexes, k_bb_indexes, k_aa_skyline, k_aa_matrix, k_ab_matrix, k_ba_matrix, k_bb_matrix,
            )
        )
    }
}
