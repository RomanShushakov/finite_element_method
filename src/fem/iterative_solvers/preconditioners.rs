use extended_matrix::{BasicOperationsTrait, FloatTrait, Position, SquareMatrix, Vector};

/// Jacobi preconditioner: z = M^{-1} r, where M = diag(A)
pub fn apply_jacobi_preconditioner<V>(
    a: &SquareMatrix<V>,
    r: &Vector<V>,
) -> Result<Vector<V>, String>
where
    V: FloatTrait<Output = V>,
{
    let a_shape = a.get_shape();
    let r_shape = r.get_shape();
    let (rows_a, cols_a) = (a_shape.0, a_shape.1);
    let (rows_r, cols_r) = (r_shape.0, r_shape.1);

    // r must be a vector and size must match A
    if cols_r != 1 && rows_r != 1 {
        return Err("Jacobi: r is not a vector".to_string());
    }

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
        // column vector
        for i in 0..n {
            let rii = *r
                .get_element_value(&Position(i, 0))
                .map_err(|e| format!("Jacobi (r col): {}", e))?;
            let aii = *a
                .get_element_value(&Position(i, i))
                .map_err(|e| format!("Jacobi (A diag): {}", e))?;

            // For SPD matrices we expect aii > 0; you may later add checks/shifts.
            z_vals.push(rii / aii);
        }
    } else {
        // row vector
        for j in 0..n {
            let rjj = *r
                .get_element_value(&Position(0, j))
                .map_err(|e| format!("Jacobi (r row): {}", e))?;
            let ajj = *a
                .get_element_value(&Position(j, j))
                .map_err(|e| format!("Jacobi (A diag): {}", e))?;
            z_vals.push(rjj / ajj);
        }
    }

    Ok(Vector::create(&z_vals))
}
