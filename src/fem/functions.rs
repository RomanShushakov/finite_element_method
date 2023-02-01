use std::collections::HashMap;

use extended_matrix::
{
    FloatTrait, Matrix, SquareMatrix, Position, BasicOperationsTrait, Vector3, VectorTrait,
};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::plate::plate4n4ip::Plate4n4ip;
use crate::fem::finite_elements::beam::beam2n1ipt::Beam2n1ipT;
use crate::fem::global_analysis::fe_stiffness::
{
    stiffness_types_number, StiffnessType, StiffnessGroupKey,
};
use crate::fem::global_analysis::fe_dof_parameter_data::global_dof;

use crate::fem::separated_matrix::SeparatedMatrix;

use crate::fem::finite_elements::functions::compare_with_tolerance;

use crate::fem::convex_hull_on_plane::{Point, convex_hull_on_plane};


pub(super) fn add_new_stiffness_sub_groups(
    stiffness_groups: &mut HashMap<StiffnessGroupKey, Vec<Position>>,
    global_group_position: usize, 
    global_group_columns_number: usize, 
    global_number_1: u32,
    global_number_2: u32
) 
    -> Result<(), String>
{
    let row = global_group_position / global_group_columns_number;
    let column = global_group_position % global_group_columns_number;

    let mut k = 0;
    while k < stiffness_types_number()
    {
        let start_row = row * global_dof();

        let row_shift_init = k / 2 * (global_dof() / 2);

        let row_shift_final = k / 2 * (global_dof() / 2) + (global_dof() / 2);

        let start_column = column * global_dof();

        let column_shift_init = k % 2 * (global_dof() / 2);

        let column_shift_final = k % 2 * (global_dof() / 2) + (global_dof() / 2);

        let mut element_positions = Vec::new();
        let mut current_row = start_row + row_shift_init;
        while current_row < start_row + row_shift_final
        {
            let mut current_column  = start_column + column_shift_init;
            while current_column < start_column + column_shift_final
            {
                element_positions.push(Position(current_row, current_column));
                current_column += 1;
            }
            current_row += 1;
        }
        let stiffness_type = StiffnessType::iterator()
            .nth(k)
            .ok_or("FEModel: Stiffness type could not be defined".to_string())?;
        let stiffness_group_key = StiffnessGroupKey { stiffness_type: *stiffness_type,
            number_1: global_number_1, number_2: global_number_2 };
        stiffness_groups.insert(stiffness_group_key, element_positions);
        k += 1;
    }
    Ok(())
}


pub(super) fn separate<V>(matrix: Matrix<V>, positions: Vec<Position>) -> Result<SeparatedMatrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let shape = matrix.get_shape();

    let k_aa_order = shape.0 - positions.len();

    let mut k_aa_elements = Vec::new();

    let mut i = 0;
    while i < shape.0
    {
        let mut j = 0;
        while j < shape.1
        {
            if positions.iter().position(|p| p.0 == i).is_none() &&
                positions.iter().position(|p| p.1 == j).is_none()
            {
                let value = matrix.get_element_value(&Position(i, j))?;
                k_aa_elements.push(*value);
            }
            j += 1;
        }
        i += 1;
    }

    let k_aa_matrix = SquareMatrix::create(k_aa_order, &k_aa_elements);

    let k_ab_rows_number = shape.0 - positions.len();

    let k_ab_columns_number = positions.len();

    let mut k_ab_elements = Vec::new();

    let mut i = 0;
    while i < shape.0
    {
        if positions.iter().position(|p| p.0 == i).is_none()
        {
            for j in 0..positions.len()
            {
                let row = i;
                let column = positions[j].1;
                if column > shape.1
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kab \
                        could not be composed!".to_string());
                }
                let value = matrix.get_element_value(&Position(row, column))?;
                k_ab_elements.push(*value);
            }
        }
        i += 1;
    }

    let k_ab_matrix = Matrix::create(k_ab_rows_number, k_ab_columns_number, &k_ab_elements);

    let k_ba_rows_number = positions.len();

    let k_ba_columns_number = shape.1 - positions.len();

    let mut k_ba_elements = Vec::new();


    for i in 0..positions.len()
    {
        let mut j = 0;
        while j < shape.1
        {
            if positions.iter().position(|p| p.1 == j).is_none()
            {
                let row = positions[i].0;
                let column = j;
                if row > shape.0
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kba \
                        could not be composed!".to_string());
                }
                let value = matrix.get_element_value(&Position(row, column))?;
                k_ba_elements.push(*value);
            }
            j += 1;
        }
    }

    let k_ba_matrix = Matrix::create(k_ba_rows_number, k_ba_columns_number, &k_ba_elements);

    let k_bb_order = positions.len();

    let mut k_bb_elements = Vec::new();

    for i in 0..positions.len()
    {
        for j in 0..positions.len()
        {
            let row = positions[i].0;
            let column = positions[j].1;
            if row > shape.0 || column > shape.1
            {
                return Err("Extended matrix: Matrix could not be separated! Matrix Kbb could \
                    not be composed!".to_string());
            }
            let value = matrix.get_element_value(&Position(row, column))?;
            k_bb_elements.push(*value);
        }
    }
    let k_bb_matrix = SquareMatrix::create(k_bb_order, &k_bb_elements);

    let separated_matrix = SeparatedMatrix::create(k_aa_matrix, k_ab_matrix, k_ba_matrix, 
        k_bb_matrix);

    Ok(separated_matrix)
}


pub fn is_points_of_quadrilateral_on_the_same_line<V>(point_1: &[V], point_2: &[V], point_3: &[V], point_4: &[V],
    tolerance: V) -> bool
    where V: FloatTrait<Output = V>,
{
    let cross_product_handle = |vector_1: &[V], vector_2: &[V]| 
        {
            [
                compare_with_tolerance(vector_1[1] * vector_2[2] - vector_1[2] * vector_2[1], tolerance),
                compare_with_tolerance(vector_1[2] * vector_2[0] - vector_1[0] * vector_2[2], tolerance),
                compare_with_tolerance(vector_1[0] * vector_2[1] - vector_1[1] * vector_2[0], tolerance),
            ]
        };

    let vector_1_2 = [point_2[0] - point_1[0], point_2[1] - point_1[1], point_2[2] - point_1[2]];
    let vector_1_4 = [point_4[0] - point_1[0], point_4[1] - point_1[1], point_4[2] - point_1[2]];

    let vector_2_1 = [point_1[0] - point_2[0], point_1[1] - point_2[1], point_1[2] - point_2[2]];
    let vector_2_3 = [point_3[0] - point_2[0], point_3[1] - point_2[1], point_3[2] - point_2[2]];

    let vector_3_2 = [point_2[0] - point_3[0], point_2[1] - point_3[1], point_2[2] - point_3[2]];
    let vector_3_4 = [point_4[0] - point_3[0], point_4[1] - point_3[1], point_4[2] - point_3[2]];

    let vector_4_3 = [point_3[0] - point_4[0], point_3[1] - point_4[1], point_3[2] - point_4[2]];
    let vector_4_1 = [point_1[0] - point_4[0], point_1[1] - point_4[1], point_1[2] - point_4[2]];

    let vector_pairs = [
        [vector_1_2, vector_1_4], [vector_2_1, vector_2_3], [vector_3_2, vector_3_4], [vector_4_3, vector_4_1],
    ];

    for vector_pair in vector_pairs
    {
        let cross_product = cross_product_handle(&vector_pair[0], &vector_pair[1]);
        if cross_product[0] == V::from(0f32) && cross_product[1] == V::from(0f32) && cross_product[2] == V::from(0f32)
        {
            return true;
        }
    }

    false
}


pub fn is_points_of_quadrilateral_on_the_same_plane<V>(point_1: &[V], point_2: &[V], point_3: &[V], point_4: &[V],
    tolerance: V) -> bool
    where V: FloatTrait<Output = V>,
{
    let cross_product_handle = |vector_1: &[V], vector_2: &[V]| 
        {
            [
                compare_with_tolerance(vector_1[1] * vector_2[2] - vector_1[2] * vector_2[1], tolerance),
                compare_with_tolerance(vector_1[2] * vector_2[0] - vector_1[0] * vector_2[2], tolerance),
                compare_with_tolerance(vector_1[0] * vector_2[1] - vector_1[1] * vector_2[0], tolerance),
            ]
        };
    let vector_3_2 = [point_2[0] - point_3[0], point_2[1] - point_3[1], point_2[2] - point_3[2]];
    let vector_3_4 = [point_4[0] - point_3[0], point_4[1] - point_3[1], point_4[2] - point_3[2]]; 

    let normal_to_plane = cross_product_handle(&vector_3_2, &vector_3_4);
    let a = normal_to_plane[0];
    let b = normal_to_plane[1];
    let c = normal_to_plane[2];
    let d = V::from(-1f32) * (a * point_3[0] + b * point_3[1] + c * point_3[2]);
    
    if compare_with_tolerance(a * point_1[0] + b * point_1[1] + c * point_1[2] + d, tolerance) != V::from(0f32)
    {
        return false
    }
    true
}


fn rotation_matrix_of_quadrilateral<V>(
    point_2: &[V], 
    point_3: &[V], 
    point_4: &[V], 
    rel_tol: V,
    abs_tol: V,
) 
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let edge_3_4_x = point_4[0] - point_3[0];
    let edge_3_4_y = point_4[1] - point_3[1];
    let edge_3_4_z = point_4[2] - point_3[2];

    let edge_3_2_x = point_2[0] - point_3[0];
    let edge_3_2_y = point_2[1] - point_3[1];
    let edge_3_2_z = point_2[2] - point_3[2];
    let normal_through_node_3_x = edge_3_4_y * edge_3_2_z - edge_3_4_z * edge_3_2_y;
    let normal_through_node_3_y = edge_3_4_z * edge_3_2_x - edge_3_4_x * edge_3_2_z;
    let normal_through_node_3_z = edge_3_4_x * edge_3_2_y - edge_3_4_y * edge_3_2_x;

    let normal_through_node_3_vector = Vector3::create(
        &[normal_through_node_3_x, normal_through_node_3_y, normal_through_node_3_z],
    );
    let normal_through_node_3_length = normal_through_node_3_vector.norm()?;
    let direction_vector = Vector3::create(
        &[V::from(0f32), V::from(0f32), normal_through_node_3_length],
    );

    let rotation_matrix = normal_through_node_3_vector
        .rotation_matrix_to_align_with_vector(&direction_vector, rel_tol, abs_tol)?;

    // let normal_through_node_3_length = ((normal_through_node_3_x).my_powi(2) + 
    //     (normal_through_node_3_y).my_powi(2) + (normal_through_node_3_z).my_powi(2)).my_sqrt();
    // let (u, v, w) = (V::from(0f32), V::from(0f32), normal_through_node_3_length);
    // let alpha = ((normal_through_node_3_x * u + normal_through_node_3_y * v + normal_through_node_3_z * w) / 
    //     (normal_through_node_3_length * normal_through_node_3_length)).my_acos();

    // let rotation_axis_coord_x = 
    //     {
    //         if normal_through_node_3_x == V::from(0f32) && 
    //             normal_through_node_3_y == V::from(0f32) && 
    //             normal_through_node_3_z != V::from(0f32) 
    //         {
    //             normal_through_node_3_length
    //         }
    //         else
    //         {
    //             normal_through_node_3_y * w - normal_through_node_3_z * v
    //         }
    //     };
    // let rotation_axis_coord_y = normal_through_node_3_z * u - normal_through_node_3_x * w;
    // let rotation_axis_coord_z = normal_through_node_3_x * v - normal_through_node_3_y * u;

    // let norm = V::from(1f32) / (rotation_axis_coord_x.my_powi(2) +
    //     rotation_axis_coord_y.my_powi(2) + rotation_axis_coord_z.my_powi(2)).my_sqrt();
    // let (x_n, y_n, z_n) = (rotation_axis_coord_x * norm,
    //     rotation_axis_coord_y * norm, rotation_axis_coord_z * norm);
    // let (c, s) = (alpha.my_cos(), alpha.my_sin());
    // let t = V::from(1f32) - c;
    // let q_11 = compare_with_tolerance(t * x_n * x_n + c, tolerance);
    // let q_12 = compare_with_tolerance(t * x_n * y_n - z_n * s, tolerance);
    // let q_13 = compare_with_tolerance(t * x_n * z_n + y_n * s, tolerance);
    // let q_21 = compare_with_tolerance(t * x_n * y_n + z_n * s, tolerance);
    // let q_22 = compare_with_tolerance(t * y_n * y_n + c, tolerance);
    // let q_23 = compare_with_tolerance(t * y_n * z_n - x_n * s, tolerance);
    // let q_31 = compare_with_tolerance(t * x_n * z_n - y_n * s, tolerance);
    // let q_32 = compare_with_tolerance(t * y_n * z_n + x_n * s, tolerance);
    // let q_33 = compare_with_tolerance(t * z_n * z_n + c, tolerance);
    
    // let rotation_matrix = SquareMatrix::create(3usize, 
    //     &[q_11, q_12, q_13, q_21, q_22, q_23, q_31, q_32, q_33]);

    Ok(rotation_matrix)
}


pub fn convex_hull_on_four_points_on_plane<V>(
    point_numbers: &[u32], 
    points: &[&[V]], 
    rel_tol: V,
    abs_tol: V,
) 
    -> Result<Vec<u32>, String>
    where V: FloatTrait<Output = V>
{
    let rotation_matrix = rotation_matrix_of_quadrilateral::<V>(
        points[1], points[2], points[3], rel_tol, abs_tol,
    )?;

    let point_1_direction_vector = vec![
        points[0][0] - points[2][0], points[0][1] - points[2][1], points[0][2] - points[2][2],
    ];
    let point_2_direction_vector = vec![
        points[1][0] - points[2][0], points[1][1] - points[2][1], points[1][2] - points[2][2],
    ];
    let point_4_direction_vector = vec![
        points[3][0] - points[2][0], points[3][1] - points[2][1], points[3][2] - points[2][2],
    ];

    let point_1_direction = Matrix::create(
        3usize, 1, &point_1_direction_vector
    );
    let point_2_direction = Matrix::create(
        3usize, 1, &point_2_direction_vector,
    );
    let point_4_direction = Matrix::create(
        3usize, 1, &point_4_direction_vector, 
    );

    let transformed_point_1_direction = rotation_matrix.multiply(&point_1_direction)?;
    let transformed_point_2_direction = rotation_matrix.multiply(&point_2_direction)?;
    let transformed_point_4_direction = rotation_matrix.multiply(&point_4_direction)?;

    let transformed_point_1_direction_x = *transformed_point_1_direction.get_element_value(&Position(0, 0))?;
    let transformed_point_1_direction_y = *transformed_point_1_direction.get_element_value(&Position(1, 0))?;
    let transformed_point_2_direction_x = *transformed_point_2_direction.get_element_value(&Position(0, 0))?;
    let transformed_point_2_direction_y = *transformed_point_2_direction.get_element_value(&Position(1, 0))?;
    let transformed_point_4_direction_x = *transformed_point_4_direction.get_element_value(&Position(0, 0))?;
    let transformed_point_4_direction_y = *transformed_point_4_direction.get_element_value(&Position(1, 0))?;

    let point_1_on_plane = Point::create(
        point_numbers[0], transformed_point_1_direction_x, transformed_point_1_direction_y,
    );
    let point_2_on_plane = Point::create(
        point_numbers[1], transformed_point_2_direction_x, transformed_point_2_direction_y,
    );
    let point_3_on_plane = Point::create(
        point_numbers[2], V::from(0f32), V::from(0f32),
    );
    let point_4_on_plane = Point::create(
        point_numbers[3], transformed_point_4_direction_x, transformed_point_4_direction_y,
    );

    let convex_hull_on_plane = convex_hull_on_plane(
        &[point_1_on_plane, point_2_on_plane, point_3_on_plane, point_4_on_plane]
    );
    
    let convex_hull_point_numbers = convex_hull_on_plane
        .iter()
        .map(|point| point.copy_number())
        .collect::<Vec<u32>>();

    Ok(convex_hull_point_numbers)
}


pub fn convert_uniformly_distributed_surface_force_to_nodal_forces<V>(
    node_1_data: (u32, V, V, V),
    node_2_data: (u32, V, V, V),
    node_3_data: (u32, V, V, V),
    node_4_data: (u32, V, V, V), 
    uniformly_distributed_surface_force_value: V,
    rel_tol: V,
    abs_tol: V,
)
    -> Result<HashMap<u32, V>, String>
    where V: FloatTrait<Output = V>
{
    let mut nodes = HashMap::new();
    nodes.insert(node_1_data.0, FENode::create(node_1_data.1, node_1_data.2, node_1_data.3));
    nodes.insert(node_2_data.0, FENode::create(node_2_data.1, node_2_data.2, node_2_data.3));
    nodes.insert(node_3_data.0, FENode::create(node_3_data.1, node_3_data.2, node_3_data.3));
    nodes.insert(node_4_data.0, FENode::create(node_4_data.1, node_4_data.2, node_4_data.3));

    let default_plate_element = Plate4n4ip::create(
        node_1_data.0, 
        node_2_data.0, 
        node_3_data.0, 
        node_4_data.0, 
        V::from(1e7f32), 
        V::from(0.3f32), 
        V::from(0.1f32), 
        V::from(0.833f32), 
        &nodes,
        rel_tol,
        abs_tol,
    )?;
    
    let nodal_forces_matrix = default_plate_element
        .convert_uniformly_distributed_surface_force_to_nodal_forces(
            uniformly_distributed_surface_force_value, &nodes, rel_tol,
    )?;

    let nodal_force_1 = *nodal_forces_matrix.get_element_value(&Position(0, 0))?;
    let nodal_force_2 = *nodal_forces_matrix.get_element_value(&Position(1, 0))?;
    let nodal_force_3 = *nodal_forces_matrix.get_element_value(&Position(2, 0))?;
    let nodal_force_4 = *nodal_forces_matrix.get_element_value(&Position(3, 0))?;

    Ok(HashMap::from([
        (node_1_data.0, nodal_force_1), (node_2_data.0, nodal_force_2), 
        (node_3_data.0, nodal_force_3), (node_4_data.0, nodal_force_4),
    ]))
}


pub fn convert_uniformly_distributed_line_force_to_nodal_forces<V>(
    node_1_data: (u32, V, V, V),
    node_2_data: (u32, V, V, V), 
    uniformly_distributed_line_force_value: V, 
    rel_tol: V,
    abs_tol: V,
) 
    -> Result<HashMap<u32, V>, String>
    where V: FloatTrait<Output = V>
{
    let mut nodes = HashMap::new();
    nodes.insert(node_1_data.0, FENode::create(node_1_data.1, node_1_data.2, node_1_data.3));
    nodes.insert(node_2_data.0, FENode::create(node_2_data.1, node_2_data.2, node_2_data.3));

    let default_beam_element = Beam2n1ipT::create(
        node_1_data.0,
        node_2_data.0,
        V::from(1e7f32),
        V::from(0.3f32),
        V::from(1f32),
        V::from(1f32),
        V::from(1f32), V::from(0f32),
        V::from(1f32),
        V::from(0.833f32), 
        [V::from(0f32), V::from(0f32), V::from(1f32)],
        &nodes,
        rel_tol,
        abs_tol,
    )?;

    let nodal_forces_matrix = default_beam_element
        .convert_uniformly_distributed_line_force_to_nodal_forces(
            uniformly_distributed_line_force_value, &nodes,
    )?;
    
    let nodal_force_1 = *nodal_forces_matrix.get_element_value(&Position(0, 0))?;
    let nodal_force_2 = *nodal_forces_matrix.get_element_value(&Position(1, 0))?;

    Ok(HashMap::from([(node_1_data.0, nodal_force_1), (node_2_data.0, nodal_force_2)]))
}
