use extended_matrix::{BasicOperationsTrait, FloatTrait, Position, SquareMatrix, Vector, VectorTrait};

/// Jacobi preconditioner: z = M^{-1} r, where M = diag(A)
pub fn apply_jacobi_preconditioner<V>(
    a: &SquareMatrix<V>,
    r: &Vector<V>,
) -> Result<Vector<V>, String>
where
    V: FloatTrait<Output = V>,
{
    r.vector_shape_conformity_check()?;

    let a_shape = a.get_shape();
    let r_shape = r.get_shape();
    let (rows_a, cols_a) = (a_shape.0, a_shape.1);
    let (rows_r, cols_r) = (r_shape.0, r_shape.1);

    let n = if cols_r == 1 { rows_r } else { cols_r };

    if rows_a != cols_a || rows_a != n {
        return Err(format!(
            "Jacobi: size mismatch: A = {:?}, r (interpreted as len {})",
            a.get_shape(),
            n
        ));
    }

    let mut z_vals = Vec::with_capacity(n);

    if cols_r == 1 {
        for i in 0..n {
            let r_ii = *r
                .get_element_value(&Position(i, 0))
                .map_err(|e| format!("Jacobi (r col): {}", e))?;
            let a_ii = *a
                .get_element_value(&Position(i, i))
                .map_err(|e| format!("Jacobi (A diag): {}", e))?;

            // For SPD matrices we expect aii > 0; you may later add checks/shifts.
            z_vals.push(r_ii / a_ii);
        }
    } else {
        for j in 0..n {
            let r_jj = *r
                .get_element_value(&Position(0, j))
                .map_err(|e| format!("Jacobi (r row): {}", e))?;
            let a_jj = *a
                .get_element_value(&Position(j, j))
                .map_err(|e| format!("Jacobi (A diag): {}", e))?;
            z_vals.push(r_jj / a_jj);
        }
    }

    Ok(Vector::create(&z_vals))
}
