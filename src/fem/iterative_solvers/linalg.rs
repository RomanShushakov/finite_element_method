use extended_matrix::{BasicOperationsTrait, FloatTrait, Position, Vector, VectorTrait};

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
    x.vector_shape_conformity_check()?;
    y.vector_shape_conformity_check()?;

    let x_shape = x.get_shape();
    let y_shape = y.get_shape();

    if x_shape != y_shape {
        return Err(format!(
            "axpy: shape mismatch: x = {:?}, y = {:?}",
            x.get_shape(),
            y.get_shape()
        ));
    }

    if x_shape.1 == 1 {
        for i in 0..x_shape.0 {
            let x_i = *x
                .get_element_value(&Position(i, 0))
                .map_err(|e| format!("axpy (x col): {}", e))?;
            let y_i = y
                .get_mut_element_value(&Position(i, 0))
                .map_err(|e| format!("axpy (y col): {}", e))?;
            *y_i = *y_i + alpha * x_i;
        }
    } else {
        for j in 0..x_shape.1 {
            let x_j = *x
                .get_element_value(&Position(0, j))
                .map_err(|e| format!("axpy (x row): {}", e))?;
            let y_j = y
                .get_mut_element_value(&Position(0, j))
                .map_err(|e| format!("axpy (y row): {}", e))?;
            *y_j = *y_j + alpha * x_j;
        }
    }

    Ok(())
}

pub fn scale<V>(x: &mut Vector<V>, alpha: V) -> Result<(), String>
where
    V: FloatTrait<Output = V>,
{
    x.vector_shape_conformity_check()?;
    let scaled_x = x.multiply_by_scalar(alpha);
    *x = scaled_x;
    Ok(())
}
