use std::collections::HashMap;

use extended_matrix::{
    BasicOperationsTrait, FloatTrait, Matrix, Position, SquareMatrix, SquareMatrixTrait, Vector,
    Vector3, VectorTrait,
};

use crate::fem::bar_2n_element_functions::find_2n_element_vector;
use crate::fem::convex_hull_on_plane::{Point, convex_hull_on_plane};
use crate::fem::math_functions::compare_with_tolerance;
use crate::fem::math_functions::{derivative_x, power_func_x};
use crate::fem::structs::Node;

pub fn is_points_of_quadrilateral_on_the_same_line<V>(
    point_1: &[V],
    point_2: &[V],
    point_3: &[V],
    point_4: &[V],
) -> bool
where
    V: FloatTrait<Output = V>,
{
    let vector_1_2 = Vector3::create(&[
        point_2[0] - point_1[0],
        point_2[1] - point_1[1],
        point_2[2] - point_1[2],
    ]);
    let vector_1_4 = Vector3::create(&[
        point_4[0] - point_1[0],
        point_4[1] - point_1[1],
        point_4[2] - point_1[2],
    ]);

    let vector_2_1 = Vector3::create(&[
        point_1[0] - point_2[0],
        point_1[1] - point_2[1],
        point_1[2] - point_2[2],
    ]);
    let vector_2_3 = Vector3::create(&[
        point_3[0] - point_2[0],
        point_3[1] - point_2[1],
        point_3[2] - point_2[2],
    ]);

    let vector_3_2 = Vector3::create(&[
        point_2[0] - point_3[0],
        point_2[1] - point_3[1],
        point_2[2] - point_3[2],
    ]);
    let vector_3_4 = Vector3::create(&[
        point_4[0] - point_3[0],
        point_4[1] - point_3[1],
        point_4[2] - point_3[2],
    ]);

    let vector_4_3 = Vector3::create(&[
        point_3[0] - point_4[0],
        point_3[1] - point_4[1],
        point_3[2] - point_4[2],
    ]);
    let vector_4_1 = Vector3::create(&[
        point_1[0] - point_4[0],
        point_1[1] - point_4[1],
        point_1[2] - point_4[2],
    ]);

    let vector_pairs = [
        [vector_1_2, vector_1_4],
        [vector_2_1, vector_2_3],
        [vector_3_2, vector_3_4],
        [vector_4_3, vector_4_1],
    ];

    for vector_pair in vector_pairs {
        if vector_pair[0]
            .cross_product(&vector_pair[1])
            .norm()
            .expect("Norm could not be found!")
            == V::from(0f32)
        {
            return true;
        }
    }

    false
}

pub fn is_points_of_quadrilateral_on_the_same_plane<V>(
    point_1: &[V],
    point_2: &[V],
    point_3: &[V],
    point_4: &[V],
    abs_tol: V,
) -> bool
where
    V: FloatTrait<Output = V>,
{
    let vector_3_2 = Vector3::create(&[
        point_2[0] - point_3[0],
        point_2[1] - point_3[1],
        point_2[2] - point_3[2],
    ]);
    let vector_3_4 = Vector3::create(&[
        point_4[0] - point_3[0],
        point_4[1] - point_3[1],
        point_4[2] - point_3[2],
    ]);

    let normal_to_plane = vector_3_2.cross_product(&vector_3_4);
    let a = *normal_to_plane
        .get_element_value(&Position(0, 0))
        .expect("Incorrect position (0, 0)!");
    let b = *normal_to_plane
        .get_element_value(&Position(1, 0))
        .expect("Incorrect position (1, 0)!");
    let c = *normal_to_plane
        .get_element_value(&Position(2, 0))
        .expect("Incorrect position (0, 0)!");
    let d = V::from(-1f32) * (a * point_3[0] + b * point_3[1] + c * point_3[2]);

    if compare_with_tolerance(
        a * point_1[0] + b * point_1[1] + c * point_1[2] + d,
        abs_tol,
    ) != V::from(0f32)
    {
        return false;
    }
    true
}

pub fn find_rotation_matrix_elements_of_quadrilateral<V>(
    point_2: &[V],
    point_3: &[V],
    point_4: &[V],
    rel_tol: V,
    abs_tol: V,
) -> Result<[V; 9], String>
where
    V: FloatTrait<Output = V>,
{
    let edge_3_4 = Vector3::create(&[
        point_4[0] - point_3[0],
        point_4[1] - point_3[1],
        point_4[2] - point_3[2],
    ]);

    let edge_3_2 = Vector3::create(&[
        point_2[0] - point_3[0],
        point_2[1] - point_3[1],
        point_2[2] - point_3[2],
    ]);

    let normal_through_node_3_vector = edge_3_4.cross_product(&edge_3_2);
    let normal_through_node_3_length = normal_through_node_3_vector.norm()?;
    let direction_vector =
        Vector3::create(&[V::from(0f32), V::from(0f32), normal_through_node_3_length]);

    let rotation_matrix = normal_through_node_3_vector.rotation_matrix_to_align_with_vector(
        &direction_vector,
        rel_tol,
        abs_tol,
    )?;

    let mut rotation_matrix_elements = [V::from(0f32); 9];
    for i in 0..9 {
        let element_value = *rotation_matrix.get_element_value(&Position(i / 3, i % 3))?;
        rotation_matrix_elements[i] = element_value;
    }

    Ok(rotation_matrix_elements)
}

pub fn convex_hull_on_four_points_on_plane<V>(
    point_numbers: &[u32],
    points: &[&[V]],
    rel_tol: V,
    abs_tol: V,
) -> Result<Vec<u32>, String>
where
    V: FloatTrait<Output = V>,
{
    let rotation_matrix_elements = find_rotation_matrix_elements_of_quadrilateral::<V>(
        points[1], points[2], points[3], rel_tol, abs_tol,
    )?;
    let rotation_matrix = Matrix::create(3, 3, &rotation_matrix_elements);

    let point_1_direction = Vector3::create(&[
        points[0][0] - points[2][0],
        points[0][1] - points[2][1],
        points[0][2] - points[2][2],
    ]);
    let point_2_direction = Vector3::create(&[
        points[1][0] - points[2][0],
        points[1][1] - points[2][1],
        points[1][2] - points[2][2],
    ]);
    let point_4_direction = Vector3::create(&[
        points[3][0] - points[2][0],
        points[3][1] - points[2][1],
        points[3][2] - points[2][2],
    ]);

    let transformed_point_1_direction = rotation_matrix.multiply(&point_1_direction)?;
    let transformed_point_2_direction = rotation_matrix.multiply(&point_2_direction)?;
    let transformed_point_4_direction = rotation_matrix.multiply(&point_4_direction)?;

    let transformed_point_1_direction_x =
        *transformed_point_1_direction.get_element_value(&Position(0, 0))?;
    let transformed_point_1_direction_y =
        *transformed_point_1_direction.get_element_value(&Position(1, 0))?;
    let transformed_point_2_direction_x =
        *transformed_point_2_direction.get_element_value(&Position(0, 0))?;
    let transformed_point_2_direction_y =
        *transformed_point_2_direction.get_element_value(&Position(1, 0))?;
    let transformed_point_4_direction_x =
        *transformed_point_4_direction.get_element_value(&Position(0, 0))?;
    let transformed_point_4_direction_y =
        *transformed_point_4_direction.get_element_value(&Position(1, 0))?;

    let point_1_on_plane = Point::create(
        point_numbers[0],
        transformed_point_1_direction_x,
        transformed_point_1_direction_y,
    );
    let point_2_on_plane = Point::create(
        point_numbers[1],
        transformed_point_2_direction_x,
        transformed_point_2_direction_y,
    );
    let point_3_on_plane = Point::create(point_numbers[2], V::from(0f32), V::from(0f32));
    let point_4_on_plane = Point::create(
        point_numbers[3],
        transformed_point_4_direction_x,
        transformed_point_4_direction_y,
    );

    let convex_hull_on_plane = convex_hull_on_plane(&[
        point_1_on_plane,
        point_2_on_plane,
        point_3_on_plane,
        point_4_on_plane,
    ]);

    let convex_hull_point_numbers = convex_hull_on_plane
        .iter()
        .map(|point| point.copy_number())
        .collect::<Vec<u32>>();

    Ok(convex_hull_point_numbers)
}

fn dx_dr<V>(x_1: V, x_2: V, x_3: V, x_4: V, r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, x_1 * V::from(0.25f32) * s, V::from(0f32), 0)
        + derivative_x(power_func_x, x_1 * V::from(0.25f32), r, 1)
        + derivative_x(power_func_x, x_1 * V::from(0.25f32) * s, r, 1)
        + derivative_x(power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, x_2 * V::from(0.25f32) * s, V::from(0f32), 0)
        - derivative_x(power_func_x, x_2 * V::from(0.25f32), r, 1)
        - derivative_x(power_func_x, x_2 * V::from(0.25f32) * s, r, 1)
        + derivative_x(power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, x_3 * V::from(0.25f32) * s, V::from(0f32), 0)
        - derivative_x(power_func_x, x_3 * V::from(0.25f32), r, 1)
        + derivative_x(power_func_x, x_3 * V::from(0.25f32) * s, r, 1)
        + derivative_x(power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, x_4 * V::from(0.25f32) * s, V::from(0f32), 0)
        + derivative_x(power_func_x, x_4 * V::from(0.25f32), r, 1)
        - derivative_x(power_func_x, x_4 * V::from(0.25f32) * s, r, 1)
}

fn dx_ds<V>(x_1: V, x_2: V, x_3: V, x_4: V, r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, x_1 * V::from(0.25f32), s, 1)
        + derivative_x(power_func_x, x_1 * V::from(0.25f32) * r, V::from(0f32), 0)
        + derivative_x(power_func_x, x_1 * V::from(0.25f32) * r, s, 1)
        + derivative_x(power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, x_2 * V::from(0.25f32), s, 1)
        - derivative_x(power_func_x, x_2 * V::from(0.25f32) * r, V::from(0f32), 0)
        - derivative_x(power_func_x, x_2 * V::from(0.25f32) * r, s, 1)
        + derivative_x(power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, x_3 * V::from(0.25f32), s, 1)
        - derivative_x(power_func_x, x_3 * V::from(0.25f32) * r, V::from(0f32), 0)
        + derivative_x(power_func_x, x_3 * V::from(0.25f32) * r, s, 1)
        + derivative_x(power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, x_4 * V::from(0.25f32), s, 1)
        + derivative_x(power_func_x, x_4 * V::from(0.25f32) * r, V::from(0f32), 0)
        - derivative_x(power_func_x, x_4 * V::from(0.25f32) * r, s, 1)
}

fn dy_dr<V>(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, y_1 * V::from(0.25f32) * s, V::from(0f32), 0)
        + derivative_x(power_func_x, y_1 * V::from(0.25f32), r, 1)
        + derivative_x(power_func_x, y_1 * V::from(0.25f32) * s, r, 1)
        + derivative_x(power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, y_2 * V::from(0.25f32) * s, V::from(0f32), 0)
        - derivative_x(power_func_x, y_2 * V::from(0.25f32), r, 1)
        - derivative_x(power_func_x, y_2 * V::from(0.25f32) * s, r, 1)
        + derivative_x(power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, y_3 * V::from(0.25f32) * s, V::from(0f32), 0)
        - derivative_x(power_func_x, y_3 * V::from(0.25f32), r, 1)
        + derivative_x(power_func_x, y_3 * V::from(0.25f32) * s, r, 1)
        + derivative_x(power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, y_4 * V::from(0.25f32) * s, V::from(0f32), 0)
        + derivative_x(power_func_x, y_4 * V::from(0.25f32), r, 1)
        - derivative_x(power_func_x, y_4 * V::from(0.25f32) * s, r, 1)
}

fn dy_ds<V>(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, y_1 * V::from(0.25f32), s, 1)
        + derivative_x(power_func_x, y_1 * V::from(0.25f32) * r, V::from(0f32), 0)
        + derivative_x(power_func_x, y_1 * V::from(0.25f32) * r, s, 1)
        + derivative_x(power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, y_2 * V::from(0.25f32), s, 1)
        - derivative_x(power_func_x, y_2 * V::from(0.25f32) * r, V::from(0f32), 0)
        - derivative_x(power_func_x, y_2 * V::from(0.25f32) * r, s, 1)
        + derivative_x(power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, y_3 * V::from(0.25f32), s, 1)
        - derivative_x(power_func_x, y_3 * V::from(0.25f32) * r, V::from(0f32), 0)
        + derivative_x(power_func_x, y_3 * V::from(0.25f32) * r, s, 1)
        + derivative_x(power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, y_4 * V::from(0.25f32), s, 1)
        + derivative_x(power_func_x, y_4 * V::from(0.25f32) * r, V::from(0f32), 0)
        - derivative_x(power_func_x, y_4 * V::from(0.25f32) * r, s, 1)
}

pub fn extract_transformed_directions_of_nodes<V>(
    node_1_number: u32,
    node_2_number: u32,
    node_3_number: u32,
    node_4_number: u32,
    nodes: &HashMap<u32, Node<V>>,
    rotation_matrix_elements: &[V; 9],
) -> Result<[[V; 3]; 3], String>
where
    V: FloatTrait<Output = V>,
{
    let rotation_matrix = SquareMatrix::create(3, rotation_matrix_elements);

    let node_1_direction = find_2n_element_vector(node_3_number, node_1_number, nodes)?;
    let node_2_direction = find_2n_element_vector(node_3_number, node_2_number, nodes)?;
    let node_4_direction = find_2n_element_vector(node_3_number, node_4_number, nodes)?;

    let transformed_node_1_direction = rotation_matrix.multiply(&node_1_direction)?;
    let transformed_node_2_direction = rotation_matrix.multiply(&node_2_direction)?;
    let transformed_node_4_direction = rotation_matrix.multiply(&node_4_direction)?;

    Ok([
        [
            *transformed_node_1_direction.get_element_value(&Position(0, 0))?,
            *transformed_node_1_direction.get_element_value(&Position(1, 0))?,
            *transformed_node_1_direction.get_element_value(&Position(2, 0))?,
        ],
        [
            *transformed_node_2_direction.get_element_value(&Position(0, 0))?,
            *transformed_node_2_direction.get_element_value(&Position(1, 0))?,
            *transformed_node_2_direction.get_element_value(&Position(2, 0))?,
        ],
        [
            *transformed_node_4_direction.get_element_value(&Position(0, 0))?,
            *transformed_node_4_direction.get_element_value(&Position(1, 0))?,
            *transformed_node_4_direction.get_element_value(&Position(2, 0))?,
        ],
    ])
}

fn jacobian_at_r_s<V>(
    node_1_number: u32,
    node_2_number: u32,
    node_3_number: u32,
    node_4_number: u32,
    r: V,
    s: V,
    nodes: &HashMap<u32, Node<V>>,
    rotation_matrix_elements: &[V; 9],
) -> Result<SquareMatrix<V>, String>
where
    V: FloatTrait<Output = V>,
{
    let transformed_directions_of_nodes = extract_transformed_directions_of_nodes(
        node_1_number,
        node_2_number,
        node_3_number,
        node_4_number,
        nodes,
        rotation_matrix_elements,
    )?;

    let transformed_node_1_direction_x = transformed_directions_of_nodes[0][0];
    let transformed_node_1_direction_y = transformed_directions_of_nodes[0][1];
    let transformed_node_2_direction_x = transformed_directions_of_nodes[1][0];
    let transformed_node_2_direction_y = transformed_directions_of_nodes[1][1];
    let transformed_node_4_direction_x = transformed_directions_of_nodes[2][0];
    let transformed_node_4_direction_y = transformed_directions_of_nodes[2][1];

    Ok(SquareMatrix::create(
        2,
        &[
            dx_dr(
                transformed_node_1_direction_x,
                transformed_node_2_direction_x,
                V::from(0f32),
                transformed_node_4_direction_x,
                r,
                s,
            ),
            dy_dr(
                transformed_node_1_direction_y,
                transformed_node_2_direction_y,
                V::from(0f32),
                transformed_node_4_direction_y,
                r,
                s,
            ),
            dx_ds(
                transformed_node_1_direction_x,
                transformed_node_2_direction_x,
                V::from(0f32),
                transformed_node_4_direction_x,
                r,
                s,
            ),
            dy_ds(
                transformed_node_1_direction_y,
                transformed_node_2_direction_y,
                V::from(0f32),
                transformed_node_4_direction_y,
                r,
                s,
            ),
        ],
    ))
}

fn inverse_jacobian_at_r_s<V>(
    node_1_number: u32,
    node_2_number: u32,
    node_3_number: u32,
    node_4_number: u32,
    r: V,
    s: V,
    nodes: &HashMap<u32, Node<V>>,
    rotation_matrix_elements: &[V; 9],
    rel_tol: V,
) -> Result<SquareMatrix<V>, String>
where
    V: FloatTrait<Output = V>,
{
    let jacobian = jacobian_at_r_s(
        node_1_number,
        node_2_number,
        node_3_number,
        node_4_number,
        r,
        s,
        nodes,
        rotation_matrix_elements,
    )?;
    let mut x = Vector::create(&vec![V::from(0f32); jacobian.get_shape().0]);
    let inverse_jacobian = jacobian.inverse(&mut x, rel_tol)?;
    Ok(inverse_jacobian)
}

pub fn determinant_of_jacobian_at_r_s<V>(
    node_1_number: u32,
    node_2_number: u32,
    node_3_number: u32,
    node_4_number: u32,
    r: V,
    s: V,
    nodes: &HashMap<u32, Node<V>>,
    rotation_matrix_elements: &[V; 9],
    rel_tol: V,
) -> Result<V, String>
where
    V: FloatTrait<Output = V>,
{
    let jacobian = jacobian_at_r_s(
        node_1_number,
        node_2_number,
        node_3_number,
        node_4_number,
        r,
        s,
        nodes,
        rotation_matrix_elements,
    )?;
    let determinant_of_jacobian = jacobian.determinant(rel_tol);
    Ok(determinant_of_jacobian)
}

fn dh1_dr<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, V::from(0.25f32) * s, V::from(0f32), 0)
        + derivative_x(power_func_x, V::from(0.25f32), r, 1)
        + derivative_x(power_func_x, V::from(0.25f32) * s, r, 1)
}

fn dh2_dr<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, V::from(0.25f32) * s, V::from(0f32), 0)
        - derivative_x(power_func_x, V::from(0.25f32), r, 1)
        - derivative_x(power_func_x, V::from(0.25f32) * s, r, 1)
}

fn dh3_dr<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, V::from(0.25f32) * s, V::from(0f32), 0)
        - derivative_x(power_func_x, V::from(0.25f32), r, 1)
        + derivative_x(power_func_x, V::from(0.25f32) * s, r, 1)
}

fn dh4_dr<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, V::from(0.25f32) * s, V::from(0f32), 0)
        + derivative_x(power_func_x, V::from(0.25f32), r, 1)
        - derivative_x(power_func_x, V::from(0.25f32) * s, r, 1)
}

fn dh1_ds<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, V::from(0.25f32), s, 1)
        + derivative_x(power_func_x, V::from(0.25f32) * r, V::from(0f32), 0)
        + derivative_x(power_func_x, V::from(0.25f32) * r, s, 1)
}

fn dh2_ds<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, V::from(0.25f32), V::from(0f32), 0)
        + derivative_x(power_func_x, V::from(0.25f32), s, 1)
        - derivative_x(power_func_x, V::from(0.25f32) * r, V::from(0f32), 0)
        - derivative_x(power_func_x, V::from(0.25f32) * r, s, 1)
}

fn dh3_ds<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, V::from(0.25f32), s, 1)
        - derivative_x(power_func_x, V::from(0.25f32) * r, V::from(0f32), 0)
        + derivative_x(power_func_x, V::from(0.25f32) * r, s, 1)
}

fn dh4_ds<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    derivative_x(power_func_x, V::from(0.25f32), V::from(0f32), 0)
        - derivative_x(power_func_x, V::from(0.25f32), s, 1)
        + derivative_x(power_func_x, V::from(0.25f32) * r, V::from(0f32), 0)
        - derivative_x(power_func_x, V::from(0.25f32) * r, s, 1)
}

pub fn h1_r_s<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    V::from(0.25f32) * (V::from(1f32) + r) * (V::from(1f32) + s)
}

pub fn h2_r_s<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    V::from(0.25f32) * (V::from(1f32) - r) * (V::from(1f32) + s)
}

pub fn h3_r_s<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    V::from(0.25f32) * (V::from(1f32) - r) * (V::from(1f32) - s)
}

pub fn h4_r_s<V>(r: V, s: V) -> V
where
    V: FloatTrait<Output = V>,
{
    V::from(0.25f32) * (V::from(1f32) + r) * (V::from(1f32) - s)
}

pub fn dh_dx_dh_dy<V>(
    node_1_number: u32,
    node_2_number: u32,
    node_3_number: u32,
    node_4_number: u32,
    r: V,
    s: V,
    nodes: &HashMap<u32, Node<V>>,
    rotation_matrix_elements: &[V; 9],
    rel_tol: V,
) -> Result<Matrix<V>, String>
where
    V: FloatTrait<Output = V>,
{
    let inverse_jacobian = inverse_jacobian_at_r_s(
        node_1_number,
        node_2_number,
        node_3_number,
        node_4_number,
        r,
        s,
        nodes,
        rotation_matrix_elements,
        rel_tol,
    )?;
    let dh_dr_dh_ds = Matrix::create(
        2,
        4,
        &[
            dh1_dr(r, s),
            dh2_dr(r, s),
            dh3_dr(r, s),
            dh4_dr(r, s),
            dh1_ds(r, s),
            dh2_ds(r, s),
            dh3_ds(r, s),
            dh4_ds(r, s),
        ],
    );
    Ok(inverse_jacobian.multiply(&dh_dr_dh_ds)?)
}
