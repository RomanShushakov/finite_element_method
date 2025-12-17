use colsol::{factorization, find_unknown};
use extended_matrix::{
    BasicOperationsTrait, CsrMatrix, FloatTrait, Matrix, Position, SquareMatrix, Vector,
};
use iterative_solvers_smpl::{pcg_block_jacobi_csr, pcg_jacobi_csr};

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

/// Builds block boundaries in local K_aa indexing for Block Jacobi.
///
/// `k_aa_indexes[i]` is the global DOF index corresponding to row/col i in K_aa.
/// Since global DOFs are laid out as `node_index * NODE_DOF + local_dof`,
/// we group consecutive rows with the same `node_index`.
pub fn build_block_starts_from_k_aa_indexes(k_aa_indexes: &[usize]) -> Vec<usize> {
    let mut starts = Vec::new();
    if k_aa_indexes.is_empty() {
        return starts;
    }

    starts.push(0);

    let mut current_node = k_aa_indexes[0] / NODE_DOF;

    for (i, &global_dof) in k_aa_indexes.iter().enumerate() {
        let node_index = global_dof / NODE_DOF;
        if node_index != current_node {
            starts.push(i); // local index where next node starts
            current_node = node_index;
        }
    }

    starts
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

    fn build_kaa_coo_from_separated_stiffness_matrix(
        &self,
        separated_stiffness_matrix: &SeparatedStiffnessMatrix<V>,
    ) -> Result<Vec<(usize, usize, V)>, String>
    where
        V: Copy,
    {
        let k_aa = separated_stiffness_matrix.get_k_aa_matrix();
        let (n, m) = (k_aa.get_shape().0, k_aa.get_shape().1);

        if n != m {
            return Err("K_aa is not square".to_string());
        }

        let mut triplets = Vec::with_capacity(k_aa.get_elements().len());

        for (pos, val) in k_aa.get_elements().iter() {
            let i = pos.0;
            let j = pos.1;
            triplets.push((i, j, *val));
        }

        if triplets.is_empty() && n > 0 {
            return Err("K_aa has no nonzero entries".to_string());
        }

        Ok(triplets)
    }

    pub fn find_ua_vector_iterative_pcg_jacobi(
        &self,
        separated_stiffness_matrix: &SeparatedStiffnessMatrix<V>,
        r_a_vector: &Vector<V>,
        u_b_vector: &Vector<V>,
        max_iter: usize,
    ) -> Result<(Vector<V>, usize), String> {
        let b_values: Vec<V> = find_b(
            r_a_vector,
            separated_stiffness_matrix.get_k_ab_matrix(),
            u_b_vector,
        )?;
        let n = b_values.len();

        let k_aa_matrix = separated_stiffness_matrix.get_k_aa_matrix();
        let csr_matrix = CsrMatrix::from_square_matrix(k_aa_matrix)
            .map_err(|e| format!("find_ua_vector_iterative: CSR conversion failed: {}", e))?;

        if csr_matrix.get_n_rows() != n {
            return Err(format!(
                "find_ua_vector_iterative: size mismatch: K_aa rows = {}, b len = {}",
                csr_matrix.get_n_rows(),
                n
            ));
        }

        let mut u_a_values = vec![V::from(0.0_f32); n];

        let iterations = pcg_jacobi_csr(
            &csr_matrix,
            &b_values,
            &mut u_a_values,
            max_iter,
            self.get_props().get_rel_tol(),
            self.get_props().get_abs_tol(),
        )
        .map_err(|e| format!("find_ua_vector_iterative: PCG failed: {}", e))?;

        Ok((Vector::create(&u_a_values), iterations))
    }

    pub fn find_ua_vector_iterative_pcg_block_jacobi(
        &self,
        separated_stiffness_matrix: &SeparatedStiffnessMatrix<V>,
        r_a_vector: &Vector<V>,
        u_b_vector: &Vector<V>,
        max_iter: usize,
    ) -> Result<(Vector<V>, usize), String> {
        // Build RHS b = r_a - K_ab * u_b (whatever your find_b does)
        let b_values: Vec<V> = find_b(
            r_a_vector,
            separated_stiffness_matrix.get_k_ab_matrix(),
            u_b_vector,
        )?;
        let n = b_values.len();

        // K_aa → CSR (local system)
        let k_aa_matrix = separated_stiffness_matrix.get_k_aa_matrix();
        let csr_matrix = CsrMatrix::from_square_matrix(k_aa_matrix).map_err(|e| {
            format!(
                "find_ua_vector_iterative_block_jacobi: CSR conversion failed: {}",
                e
            )
        })?;

        if csr_matrix.get_n_rows() != n {
            return Err(format!(
                "find_ua_vector_iterative_block_jacobi: size mismatch: K_aa rows = {}, b len = {}",
                csr_matrix.get_n_rows(),
                n
            ));
        }

        // Build block boundaries (local indices in K_aa)
        let block_starts =
            build_block_starts_from_k_aa_indexes(separated_stiffness_matrix.get_k_aa_indexes());

        // Solve
        let mut u_a_values = vec![V::from(0.0_f32); n];
        let iterations = pcg_block_jacobi_csr(
            &csr_matrix,
            &b_values,
            &mut u_a_values,
            max_iter,
            self.get_props().get_rel_tol(),
            self.get_props().get_abs_tol(),
            &block_starts,
        )
        .map_err(|e| format!("find_ua_vector_iterative_block_jacobi: PCG failed: {}", e))?;

        Ok((Vector::create(&u_a_values), iterations))
    }

    pub fn find_ua_vector_iterative_pcg_block_jacobi_sparse(
        &self,
        separated_stiffness_matrix: &SeparatedStiffnessMatrix<V>,
        r_a_vector: &Vector<V>,
        u_b_vector: &Vector<V>,
        max_iter: usize,
    ) -> Result<(Vector<V>, usize), String> {
        let b_values: Vec<V> = find_b(
            r_a_vector,
            separated_stiffness_matrix.get_k_ab_matrix(),
            u_b_vector,
        )?;
        let n = b_values.len();

        let triplets =
            self.build_kaa_coo_from_separated_stiffness_matrix(separated_stiffness_matrix)?;
        let csr = CsrMatrix::from_coo(n, n, &triplets)
            .map_err(|e| format!("CSR from COO failed: {}", e))?;

        let block_starts =
            build_block_starts_from_k_aa_indexes(separated_stiffness_matrix.get_k_aa_indexes());
        let mut u_a_values = vec![V::from(0.0_f32); n];

        let iterations = pcg_block_jacobi_csr(
            &csr,
            &b_values,
            &mut u_a_values,
            max_iter,
            self.get_props().get_rel_tol(),
            self.get_props().get_abs_tol(),
            &block_starts,
        )
        .map_err(|e| format!("PCG(BlockJacobi) failed: {}", e))?;

        Ok((Vector::create(&u_a_values), iterations))
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
