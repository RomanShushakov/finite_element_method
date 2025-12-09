use extended_matrix::{FloatTrait, Vector, VectorTrait, BasicOperationsTrait, Position};

pub(crate) fn dot<V>(x: &Vector<V>, y: &Vector<V>) -> Result<V, String>
where
    V: FloatTrait<Output = V>,
{
    x.dot_product(y)
}

pub(crate) fn norm2<V>(x: &Vector<V>) -> Result<V, String>
where
    V: FloatTrait<Output = V>,
{
    x.norm()
}

pub fn axpy<V>(y: &mut Vector<V>, alpha: V, x: &Vector<V>) -> Result<(), String>
where
    V: FloatTrait<Output = V>,
{
    let x_shape = y.get_shape();
    let y_shape = x.get_shape();

    if x_shape != y_shape {
        return Err(format!(
            "axpy: shape mismatch: x = {:?}, y = {:?}",
            x.get_shape(),
            y.get_shape()
        ));
    }

    // Treat both row and column vectors
    if x_shape.1 == 1 {
        // column vector
        for i in 0..x_shape.0 {
            let xi = *x
                .get_element_value(&Position(i, 0))
                .map_err(|e| format!("axpy (x col): {}", e))?;
            let yi = y
                .get_mut_element_value(&Position(i, 0))
                .map_err(|e| format!("axpy (y col): {}", e))?;
            *yi = *yi + alpha * xi;
        }
    } else if x_shape.0 == 1 {
        // row vector
        for j in 0..x_shape.1 {
            let xi = *x
                .get_element_value(&Position(0, j))
                .map_err(|e| format!("axpy (x row): {}", e))?;
            let yi = y
                .get_mut_element_value(&Position(0, j))
                .map_err(|e| format!("axpy (y row): {}", e))?;
            *yi = *yi + alpha * xi;
        }
    } else {
        return Err("axpy: not a vector (neither row nor column)".to_string());
    }

    Ok(())
}


pub fn scale<V>(x: &mut Vector<V>, alpha: V) -> Result<(), String>
where
    V: FloatTrait<Output = V>,
{
    let x_shape = x.get_shape();

    if x_shape.1 == 1 {
        for i in 0..x_shape.0 {
            let xi = x
                .get_mut_element_value(&Position(i, 0))
                .map_err(|e| format!("scale (col): {}", e))?;
            *xi = *xi * alpha;
        }
    } else if x_shape.0 == 1 {
        for j in 0..x_shape.1 {
            let xi = x
                .get_mut_element_value(&Position(0, j))
                .map_err(|e| format!("scale (row): {}", e))?;
            *xi = *xi * alpha;
        }
    } else {
        return Err("scale: not a vector".to_string());
    }

    Ok(())
}
