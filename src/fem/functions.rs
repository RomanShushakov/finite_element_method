use std::ops::{Div, Rem, Mul, Add, AddAssign, Sub, SubAssign, MulAssign};
use std::fmt::Debug;
use std::hash::Hash;
use std::collections::HashMap;

use extended_matrix::matrix_element_position::MatrixElementPosition;
use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::{conversion_uint_into_usize, matrix_element_value_extractor};

use extended_matrix_float::MyFloatTrait;

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::plate::plate4n4ip::Plate4n4ip;
use crate::fem::finite_elements::beam::beam2n1ipt::Beam2n1ipT;
use crate::fem::global_analysis::fe_stiffness::
{
    stiffness_types_number, StiffnessGroup, StiffnessType, StiffnessGroupKey
};
use crate::fem::global_analysis::fe_dof_parameter_data::global_dof;

use crate::fem::separated_matrix::SeparatedMatrix;

use crate::fem::finite_elements::functions::compare_with_tolerance;

use crate::fem::convex_hull_on_plane::{Point, convex_hull_on_plane};


pub(super) fn add_new_stiffness_sub_groups<'a, T>(
    stiffness_groups: &mut HashMap<StiffnessGroupKey<T>, Vec<MatrixElementPosition<T>>>,
    global_group_position: T, global_group_columns_number: T, global_number_1: T,
    global_number_2: T) -> Result<(), &'a str>
    where T: Copy + Debug + Div<Output = T> + Rem<Output = T> + Mul<Output = T> + Add<Output = T> +
             PartialOrd + AddAssign + From<u8> + Eq + Hash + SubAssign
{
    let row = global_group_position / global_group_columns_number;
    let column = global_group_position % global_group_columns_number;

    let mut k = T::from(0u8);
    while k < stiffness_types_number()
    {
        let start_row = row * global_dof::<T>();

        let row_shift_init = k / T::from(2u8) * (global_dof::<T>() / T::from(2u8));

        let row_shift_final = k / T::from(2u8) * (global_dof::<T>() / T::from(2u8)) +
            (global_dof::<T>() / T::from(2u8));

        let start_column = column * global_dof::<T>();

        let column_shift_init = k % T::from(2u8) * (global_dof::<T>() / T::from(2u8));

        let column_shift_final = k % T::from(2u8) * (global_dof::<T>() / T::from(2u8)) +
            (global_dof::<T>() / T::from(2u8));

        let mut element_positions = Vec::new();
        let mut current_row = start_row + row_shift_init;
        while current_row < start_row + row_shift_final
        {
            let mut current_column  = start_column + column_shift_init;
            while current_column < start_column + column_shift_final
            {
                element_positions.push(
                    MatrixElementPosition::create(current_row, current_column));
                current_column += T::from(1u8);
            }
            current_row += T::from(1u8);
        }
        let converted_index = conversion_uint_into_usize(k);
        let stiffness_type = StiffnessType::iterator()
            .nth(converted_index)
            .ok_or("FEModel: Stiffness type could not be defined")?;
        let stiffness_group_key = StiffnessGroupKey { stiffness_type: *stiffness_type,
            number_1: global_number_1, number_2: global_number_2 };
        stiffness_groups.insert(stiffness_group_key, element_positions);
        k += T::from(1u8);
    }
    Ok(())
}


pub(super) fn separate<T, V>(matrix: ExtendedMatrix<T, V>, positions: Vec<MatrixElementPosition<T>>,
    tolerance: V) -> Result<SeparatedMatrix<T, V>, String>
    where T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Copy + Debug + Eq + Hash + SubAssign + PartialOrd + AddAssign +
             From<u8> + Ord + 'static,
          V: Add<Output = V> + Mul<Output = V> + Sub<Output = V> + Div<Output = V> + Copy + Debug +
             PartialEq + AddAssign + MulAssign + SubAssign + Into<f64> + From<f32> + PartialOrd +
             MyFloatTrait + 'static
{
    let shape = matrix.copy_shape();

    let mut converted_positions_length = T::from(0u8);
    (0..positions.len()).for_each(|_| converted_positions_length += T::from(1u8));

    let k_aa_rows_number = shape.0 - converted_positions_length;

    let k_aa_columns_number = shape.1 - converted_positions_length;

    let mut k_aa_elements = Vec::new();

    let mut i = T::from(0u8);
    while i < shape.0
    {
        let mut j = T::from(0u8);
        while j < shape.1
        {
            if positions.iter().position(|p| *p.ref_row() == i).is_none() &&
                positions.iter().position(|p| *p.ref_column() == j).is_none()
            {
                let value = matrix_element_value_extractor(i, j, &matrix)?;
                k_aa_elements.push(value);
            }
            j += T::from(1u8);
        }
        i += T::from(1u8);
    }

    let k_aa_matrix = ExtendedMatrix::create(k_aa_rows_number,
        k_aa_columns_number, k_aa_elements, tolerance)?;

    let k_ab_rows_number = shape.0 - converted_positions_length;

    let k_ab_columns_number = converted_positions_length;

    let mut k_ab_elements = Vec::new();

    let mut i = T::from(0u8);
    while i < shape.0
    {
        if positions.iter().position(|p| *p.ref_row() == i).is_none()
        {
            for j in 0..positions.len()
            {
                let row = i;
                let column = positions[j].ref_column();
                if *column > shape.1
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kab \
                        could not be composed!".to_string());
                }
                let value = matrix_element_value_extractor(row, *column, &matrix)?;
                k_ab_elements.push(value);
            }
        }
        i += T::from(1u8);
    }

    let k_ab_matrix = ExtendedMatrix::create(k_ab_rows_number,
        k_ab_columns_number, k_ab_elements, tolerance)?;

    let k_ba_rows_number = converted_positions_length;

    let k_ba_columns_number = shape.1 - converted_positions_length;

    let mut k_ba_elements = Vec::new();


    for i in 0..positions.len()
    {
        let mut j = T::from(0u8);
        while j < shape.1
        {
            if positions.iter().position(|p| *p.ref_column() == j).is_none()
            {
                let row = positions[i].ref_row();
                let column = j;
                if *row > shape.0
                {
                    return Err("Extended matrix: Matrix could not be separated! Matrix Kba \
                        could not be composed!".to_string());
                }
                let value = matrix_element_value_extractor(*row, column, &matrix)?;
                k_ba_elements.push(value);
            }
            j += T::from(1u8);
        }
    }

    let k_ba_matrix = ExtendedMatrix::create(k_ba_rows_number,
        k_ba_columns_number, k_ba_elements, tolerance)?;

    let k_bb_rows_number = converted_positions_length;

    let k_bb_columns_number = converted_positions_length;

    let mut k_bb_elements = Vec::new();

    for i in 0..positions.len()
    {
        for j in 0..positions.len()
        {
            let row = positions[i].ref_row();
            let column = positions[j].ref_column();
            if *row > shape.0 || *column > shape.1
            {
                return Err("Extended matrix: Matrix could not be separated! Matrix Kbb could \
                    not be composed!".to_string());
            }
            let value = matrix_element_value_extractor(*row, *column, &matrix)?;
            k_bb_elements.push(value);
        }
    }
    let k_bb_matrix = ExtendedMatrix::create(k_bb_rows_number,
        k_bb_columns_number, k_bb_elements, tolerance)?;

    let separated_matrix = SeparatedMatrix::create(k_aa_matrix,
        k_ab_matrix, k_ba_matrix, k_bb_matrix);

    Ok(separated_matrix)
}


pub fn is_points_of_quadrilateral_on_the_same_line<V>(point_1: &[V], point_2: &[V], point_3: &[V], point_4: &[V],
    tolerance: V) -> bool
    where V: Copy + Add<Output = V> + Sub<Output = V> + Mul<Output = V> + From<f32> + MyFloatTrait + PartialOrd,
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
    where V: Copy + Debug + Add<Output = V> + Sub<Output = V> + Mul<Output = V> + PartialOrd + MyFloatTrait + 
             From<f32>,
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



fn rotation_matrix_of_quadrilateral<T, V>(point_2: &[V], point_3: &[V], point_4: &[V], 
    tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    where T: Copy + Debug + From<u8> + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> + 
             Rem<Output = T> + Hash + AddAssign + SubAssign + Ord + 'static,
          V: Copy + Debug + Add<Output = V> + Sub<Output = V> + Mul<Output = V> + PartialOrd + MyFloatTrait + 
             From<f32> + Div<Output = V> + AddAssign + SubAssign + MulAssign + Into<f64> + 'static,
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

    let normal_through_node_3_length = ((normal_through_node_3_x).my_powi(2) + 
        (normal_through_node_3_y).my_powi(2) + (normal_through_node_3_z).my_powi(2)).my_sqrt();
    let (u, v, w) = (V::from(0f32), V::from(0f32), normal_through_node_3_length);
    let alpha = ((normal_through_node_3_x * u + normal_through_node_3_y * v + normal_through_node_3_z * w) / 
        (normal_through_node_3_length * normal_through_node_3_length)).my_acos();

    let rotation_axis_coord_x = 
        {
            if normal_through_node_3_x == V::from(0f32) && 
                normal_through_node_3_y == V::from(0f32) && 
                normal_through_node_3_z != V::from(0f32) 
            {
                normal_through_node_3_length
            }
            else
            {
                normal_through_node_3_y * w - normal_through_node_3_z * v
            }
        };
    let rotation_axis_coord_y = normal_through_node_3_z * u - normal_through_node_3_x * w;
    let rotation_axis_coord_z = normal_through_node_3_x * v - normal_through_node_3_y * u;

    let norm = V::from(1f32) / (rotation_axis_coord_x.my_powi(2) +
        rotation_axis_coord_y.my_powi(2) + rotation_axis_coord_z.my_powi(2)).my_sqrt();
    let (x_n, y_n, z_n) = (rotation_axis_coord_x * norm,
        rotation_axis_coord_y * norm, rotation_axis_coord_z * norm);
    let (c, s) = (alpha.my_cos(), alpha.my_sin());
    let t = V::from(1f32) - c;
    let q_11 = compare_with_tolerance(t * x_n * x_n + c, tolerance);
    let q_12 = compare_with_tolerance(t * x_n * y_n - z_n * s, tolerance);
    let q_13 = compare_with_tolerance(t * x_n * z_n + y_n * s, tolerance);
    let q_21 = compare_with_tolerance(t * x_n * y_n + z_n * s, tolerance);
    let q_22 = compare_with_tolerance(t * y_n * y_n + c, tolerance);
    let q_23 = compare_with_tolerance(t * y_n * z_n - x_n * s, tolerance);
    let q_31 = compare_with_tolerance(t * x_n * z_n - y_n * s, tolerance);
    let q_32 = compare_with_tolerance(t * y_n * z_n + x_n * s, tolerance);
    let q_33 = compare_with_tolerance(t * z_n * z_n + c, tolerance);
    
    let rotation_matrix = ExtendedMatrix::create(T::from(3u8), T::from(3u8),
        vec![q_11, q_12, q_13, q_21, q_22, q_23, q_31, q_32, q_33], tolerance)?;
    Ok(rotation_matrix)
}


pub fn convex_hull_on_four_points_on_plane<T, V>(point_numbers: &[T], points: &[&[V]], tolerance: V) 
    -> Result<Vec<T>, String>
    where T: Debug + Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> + Rem<Output = T> + 
             AddAssign + SubAssign + From<u8> + Hash + Ord + 'static,
          V: Debug + Copy + Add<Output = V> + Sub<Output = V> + Mul<Output = V> + Div<Output = V> + PartialOrd + 
             MyFloatTrait + From<f32> + Into<f64> + AddAssign + SubAssign + MulAssign + 'static,
{
    let rotation_matrix = rotation_matrix_of_quadrilateral::<T, V>(points[1], points[2], 
        points[3], tolerance)?;

    let point_1_direction_vector = vec![
        points[0][0] - points[2][0], points[0][1] - points[2][1], points[0][2] - points[2][2],
    ];
    let point_2_direction_vector = vec![
        points[1][0] - points[2][0], points[1][1] - points[2][1], points[1][2] - points[2][2],
    ];
    let point_4_direction_vector = vec![
        points[3][0] - points[2][0], points[3][1] - points[2][1], points[3][2] - points[2][2],
    ];

    let point_1_direction = ExtendedMatrix::create(T::from(3u8),
        T::from(1u8), point_1_direction_vector, tolerance)?;
    let point_2_direction = ExtendedMatrix::create(T::from(3u8),
        T::from(1u8), point_2_direction_vector, tolerance)?;
    let point_4_direction = ExtendedMatrix::create(T::from(3u8),
        T::from(1u8), point_4_direction_vector, tolerance)?;

    let transformed_point_1_direction = rotation_matrix.multiply_by_matrix(&point_1_direction)?;
    let transformed_point_2_direction = rotation_matrix.multiply_by_matrix(&point_2_direction)?;
    let transformed_point_4_direction = rotation_matrix.multiply_by_matrix(&point_4_direction)?;

    let transformed_point_1_direction_x =
        matrix_element_value_extractor(T::from(0u8), T::from(0u8), &transformed_point_1_direction)?;
    let transformed_point_1_direction_y =
        matrix_element_value_extractor(T::from(1u8), T::from(0u8), &transformed_point_1_direction)?;
    let transformed_point_2_direction_x =
        matrix_element_value_extractor(T::from(0u8), T::from(0u8), &transformed_point_2_direction)?;
    let transformed_point_2_direction_y =
        matrix_element_value_extractor(T::from(1u8), T::from(0u8), &transformed_point_2_direction)?;
    let transformed_point_4_direction_x =
        matrix_element_value_extractor(T::from(0u8), T::from(0u8), &transformed_point_4_direction)?;
    let transformed_point_4_direction_y =
        matrix_element_value_extractor(T::from(1u8), T::from(0u8), &transformed_point_4_direction)?;

    let point_1_on_plane = Point::create(point_numbers[0], transformed_point_1_direction_x, 
        transformed_point_1_direction_y);
    let point_2_on_plane = Point::create(point_numbers[1], transformed_point_2_direction_x, 
        transformed_point_2_direction_y);
    let point_3_on_plane = Point::create(point_numbers[2], V::from(0f32), V::from(0f32));
    let point_4_on_plane = Point::create(point_numbers[3], transformed_point_4_direction_x, 
        transformed_point_4_direction_y);

    let convex_hull_on_plane = convex_hull_on_plane(
        &[point_1_on_plane, point_2_on_plane, point_3_on_plane, point_4_on_plane]);
    
    let convex_hull_point_numbers = convex_hull_on_plane.iter().map(|point| point.copy_number())
        .collect::<Vec<T>>();

    Ok(convex_hull_point_numbers)
}


pub fn convert_uniformly_distributed_surface_force_to_nodal_forces<T, V>(node_1_data: (T, V, V, V),
    node_2_data: (T, V, V, V), node_3_data: (T, V, V, V), node_4_data: (T, V, V, V), 
    uniformly_distributed_surface_force_value: V, tolerance: V) -> Result<HashMap<T, V>, String>
    where T: Debug + Copy + Hash + SubAssign + Mul<Output = T> + AddAssign + From<u8> + Ord + 
             Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + 'static,
          V: Debug + Copy + PartialEq + Sub<Output = V> + Mul<Output = V> + From<f32> + Add<Output = V> +
             Div<Output = V> + AddAssign + MulAssign + SubAssign + MyFloatTrait + PartialOrd + Into<f64> + 
             'static,
{
    let mut nodes = HashMap::new();
    nodes.insert(node_1_data.0, FENode::create(node_1_data.1, node_1_data.2, node_1_data.3));
    nodes.insert(node_2_data.0, FENode::create(node_2_data.1, node_2_data.2, node_2_data.3));
    nodes.insert(node_3_data.0, FENode::create(node_3_data.1, node_3_data.2, node_3_data.3));
    nodes.insert(node_4_data.0, FENode::create(node_4_data.1, node_4_data.2, node_4_data.3));

    let default_plate_element = Plate4n4ip::create(node_1_data.0, node_2_data.0, 
        node_3_data.0, node_4_data.0, V::from(1e7f32), V::from(0.3f32), 
        V::from(0.1f32), V::from(0.833f32), tolerance, &nodes)?;
    
    let nodal_forces_matrix = 
        default_plate_element.convert_uniformly_distributed_surface_force_to_nodal_forces(
            uniformly_distributed_surface_force_value, &nodes, tolerance)?;

    let nodal_force_1 = matrix_element_value_extractor(T::from(0u8), T::from(0u8), &nodal_forces_matrix)?;
    let nodal_force_2 = matrix_element_value_extractor(T::from(1u8), T::from(0u8), &nodal_forces_matrix)?;
    let nodal_force_3 = matrix_element_value_extractor(T::from(2u8), T::from(0u8), &nodal_forces_matrix)?;
    let nodal_force_4 = matrix_element_value_extractor(T::from(3u8), T::from(0u8), &nodal_forces_matrix)?;

    Ok(HashMap::from([
        (node_1_data.0, nodal_force_1), (node_2_data.0, nodal_force_2), 
        (node_3_data.0, nodal_force_3), (node_4_data.0, nodal_force_4),
    ]))
}


pub fn convert_uniformly_distributed_line_force_to_nodal_forces<T, V>(node_1_data: (T, V, V, V),
    node_2_data: (T, V, V, V), uniformly_distributed_line_force_value: V, tolerance: V) 
    -> Result<HashMap<T, V>, String>
    where T: Debug + Copy + Hash + SubAssign + Mul<Output = T> + AddAssign + From<u8> + Ord + 
             Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + 'static,
          V: Debug + Copy + PartialEq + Sub<Output = V> + Mul<Output = V> + From<f32> + Add<Output = V> +
             Div<Output = V> + AddAssign + MulAssign + SubAssign + MyFloatTrait + PartialOrd + Into<f64> + 
             MyFloatTrait<Other = V> + 'static,
{
    let mut nodes = HashMap::new();
    nodes.insert(node_1_data.0, FENode::create(node_1_data.1, node_1_data.2, node_1_data.3));
    nodes.insert(node_2_data.0, FENode::create(node_2_data.1, node_2_data.2, node_2_data.3));

    let default_beam_element = Beam2n1ipT::create(node_1_data.0, 
        node_2_data.0, V::from(1e7f32), V::from(0.3f32), V::from(1f32),
        V::from(1f32), V::from(1f32), V::from(0f32), V::from(1f32), V::from(0.833f32), 
        [V::from(0f32), V::from(0f32), V::from(1f32)], tolerance, &nodes)?;

    let nodal_forces_matrix = 
        default_beam_element.convert_uniformly_distributed_line_force_to_nodal_forces(
            uniformly_distributed_line_force_value, &nodes, tolerance)?;
    
    let nodal_force_1 = matrix_element_value_extractor(T::from(0u8), T::from(0u8), &nodal_forces_matrix)?;
    let nodal_force_2 = matrix_element_value_extractor(T::from(1u8), T::from(0u8), &nodal_forces_matrix)?;

    Ok(HashMap::from([(node_1_data.0, nodal_force_1), (node_2_data.0, nodal_force_2)]))
}
