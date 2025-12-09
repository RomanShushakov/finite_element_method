use colsol::{factorization, find_unknown};
use extended_matrix::{BasicOperationsTrait, FloatTrait, Matrix, Position, SquareMatrix, Vector};

use crate::DOFParameter;
use crate::fem::FEM;
use crate::fem::structs::{NODE_DOF, SeparatedStiffnessMatrix};

fn find_b<V>(
    r_a_vector: &Vector<V>,
    k_ab_matrix: &Matrix<V>,
    u_b_vector: &Vector<V>,
) -> Result<Vec<V>, String>
where
    V: FloatTrait<Output = V>,
{
    let mut b = Vec::new();
    let b_vector = r_a_vector.subtract(&k_ab_matrix.multiply(u_b_vector)?)?;
    for i in 0..b_vector.get_shape().0 {
        b.push(*b_vector.get_element_value(&Position(i, 0))?);
    }

    Ok(b)
}

fn convert_k_aa_into_compacted_form<V>(
    k_aa_matrix: &SquareMatrix<V>,
    k_aa_skyline: &Vec<usize>,
) -> Result<(Vec<V>, Vec<i64>), String>
where
    V: FloatTrait<Output = V>,
{
    let mut a = Vec::new();
    let mut maxa = Vec::new();

    let mut count = 0;
    for i in 0..k_aa_skyline.len() {
        a.push(*k_aa_matrix.get_element_value(&Position(i, i))?);
        maxa.push(count);
        count += 1;

        let m = k_aa_skyline[i];
        if m != 0 {
            for j in 1..m + 1 {
                a.push(*k_aa_matrix.get_element_value(&Position(i - j, i))?);
                count += 1;
            }
        }
    }
    maxa.push(count);

    Ok((a, maxa))
}

fn find_r_r<V>(
    k_ba_matrix: &Matrix<V>,
    u_a_vector: &Vector<V>,
    k_bb_matrix: &SquareMatrix<V>,
    u_b_vector: &Vector<V>,
    r_b_vector: Vector<V>,
) -> Result<Vec<V>, String>
where
    V: FloatTrait<Output = V>,
{
    let mut r_r = Vec::new();
    let r_r_vector = k_ba_matrix
        .multiply(u_a_vector)?
        .add(&k_bb_matrix.multiply(u_b_vector)?)?
        .subtract(&r_b_vector)?;
    for i in 0..r_r_vector.get_shape().0 {
        r_r.push(*r_r_vector.get_element_value(&Position(i, 0))?);
    }
    Ok(r_r)
}

impl<V> FEM<V>
where
    V: FloatTrait<Output = V>,
{
    pub fn find_ua_vector(
        &self,
        separated_stiffness_matrix: &SeparatedStiffnessMatrix<V>,
        r_a_vector: &Vector<V>,
        u_b_vector: &Vector<V>,
    ) -> Result<Vector<V>, String> {
        let (mut a, maxa) = convert_k_aa_into_compacted_form(
            separated_stiffness_matrix.get_k_aa_matrix(),
            separated_stiffness_matrix.get_k_aa_skyline(),
        )?;

        let mut b = find_b(
            &r_a_vector,
            separated_stiffness_matrix.get_k_ab_matrix(),
            &u_b_vector,
        )?;

        let nn = b.len() as i64;
        factorization(&mut a, nn, &maxa)?;
        find_unknown(&a, &mut b, nn, &maxa);

        Ok(Vector::create(&b))
    }

    pub fn find_r_r_vector(
        &self,
        separated_stiffness_matrix: &SeparatedStiffnessMatrix<V>,
        u_a_vector: &Vector<V>,
        u_b_vector: &Vector<V>,
    ) -> Result<Vector<V>, String> {
        let mut r_b_vector =
            Vector::create(&vec![
                V::from(0f32);
                separated_stiffness_matrix.get_k_bb_indexes().len()
            ]);
        for i in 0..separated_stiffness_matrix.get_k_bb_indexes().len() {
            *r_b_vector.get_mut_element_value(&Position(i, 0))? =
                *self.get_forces_vector().get_element_value(&Position(
                    separated_stiffness_matrix.get_k_bb_indexes()[i],
                    0,
                ))?;
        }
        let r_r = find_r_r(
            separated_stiffness_matrix.get_k_ba_matrix(),
            u_a_vector,
            separated_stiffness_matrix.get_k_bb_matrix(),
            u_b_vector,
            r_b_vector,
        )?;

        Ok(Vector::create(&r_r))
    }

    pub fn compose_global_analysis_result(
        &mut self,
        k_aa_indexes: &[usize],
        k_bb_indexes: &[usize],
        u_a_vector: &Vector<V>,
        r_r_vector: &Vector<V>,
    ) -> Result<(), String> {
        for i in 0..k_aa_indexes.len() {
            let index = k_aa_indexes[i];
            *self
                .get_mut_displacements_vector()
                .get_mut_element_value(&Position(index, 0))? =
                *u_a_vector.get_element_value(&Position(i, 0))?
        }

        for i in 0..k_bb_indexes.len() {
            let index = k_bb_indexes[i];
            *self
                .get_mut_forces_vector()
                .get_mut_element_value(&Position(index, 0))? =
                *r_r_vector.get_element_value(&Position(i, 0))?
        }

        Ok(())
    }

    pub fn extract_global_analysis_result(&self) -> Result<Vec<(u32, DOFParameter, V, V)>, String> {
        let mut global_analysis_result = Vec::new();

        for (node_number, node) in self.get_nodes() {
            let index = node.get_index() * NODE_DOF;
            for i in 0..NODE_DOF {
                let dof_parameter = DOFParameter::from_usize(i);
                let displacement_value = self
                    .get_displacements_vector()
                    .get_element_value(&Position(index + i, 0))?;
                let load_value = self
                    .get_forces_vector()
                    .get_element_value(&Position(index + i, 0))?;
                global_analysis_result.push((
                    *node_number,
                    dof_parameter,
                    *displacement_value,
                    *load_value,
                ));
            }
        }

        Ok(global_analysis_result)
    }
}
