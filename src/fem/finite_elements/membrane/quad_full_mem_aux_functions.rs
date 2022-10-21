use std::collections::HashMap;

use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::matrix_element_value_extractor;
use extended_matrix::traits::{UIntTrait, FloatTrait};

use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData, GlobalDOFParameter};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::membrane::consts::{MEM4N4IP_NODES_NUMBER, MEMBRANE_NODE_DOF};
use crate::fem::finite_elements::functions::
{
    compare_with_tolerance, extract_unique_elements_of_rotation_matrix
};


pub struct QuadFullMemAuxFunctions<T, V>(T, V);


impl<T, V> QuadFullMemAuxFunctions<T, V>
    where T: UIntTrait<Output = T>,
          V: FloatTrait<Output = V, Other = V>
{
    pub fn nodes_number() -> T
    {
        let mut n = T::from(0u8);
        (0..MEM4N4IP_NODES_NUMBER).for_each(|_| n += T::from(1u8));
        n
    }


    pub fn node_dof() -> T
    {
        let mut n = T::from(0u8);
        (0..MEMBRANE_NODE_DOF).for_each(|_| n += T::from(1u8));
        n
    }


    pub fn rotation_matrix(node_2_number: T, node_3_number: T, node_4_number: T, 
        tolerance: V, ref_nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let node_2 = ref_nodes.get(&node_2_number).unwrap();
        let node_3 = ref_nodes.get(&node_3_number).unwrap();
        let node_4 = ref_nodes.get(&node_4_number).unwrap();

        let edge_3_4_x = node_4.copy_x() - node_3.copy_x();
        let edge_3_4_y = node_4.copy_y() - node_3.copy_y();
        let edge_3_4_z = node_4.copy_z() - node_3.copy_z();

        let edge_3_2_x = node_2.copy_x() - node_3.copy_x();
        let edge_3_2_y = node_2.copy_y() - node_3.copy_y();
        let edge_3_2_z = node_2.copy_z() - node_3.copy_z();
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
        
        let rotation_matrix = ExtendedMatrix::create(
            QuadFullMemAuxFunctions::<T, V>::nodes_number() * QuadFullMemAuxFunctions::<T, V>::node_dof(),
            QuadFullMemAuxFunctions::<T, V>::nodes_number() * QuadFullMemAuxFunctions::<T, V>::node_dof(),
            vec![
                [q_11, q_12, q_13], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [q_21, q_22, q_23], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [q_31, q_32, q_33], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],

                [V::from(0f32); MEMBRANE_NODE_DOF], [q_11, q_12, q_13],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [q_21, q_22, q_23],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [q_31, q_32, q_33],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],

                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [q_11, q_12, q_13], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [q_21, q_22, q_23], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [q_31, q_32, q_33], [V::from(0f32); MEMBRANE_NODE_DOF],

                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [q_11, q_12, q_13],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [q_21, q_22, q_23],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [q_31, q_32, q_33],
            ].concat(),
            tolerance,
        )?;
        Ok(rotation_matrix)
    }


    fn power_func_x(a: V, x: V, n: i32) -> V
    {
        (0..n).fold(a, |acc, _| acc * x)
    }


    fn derivative_x(f: fn(V, V, i32) -> V, a: V, x: V, n: i32) -> V
    {
        let mut converted_n = V::from(0f32);
        (0..n).for_each(|_| converted_n += V::from(1f32));

        f(a * converted_n, x, n - 1)
    }


    fn dx_dr(x_1: V, x_2: V, x_3: V, x_4: V, r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), r, 1) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), r, 1) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * s, r, 1)
    }


    fn dx_ds(x_1: V, x_2: V, x_3: V, x_4: V, r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), s, 1) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), s, 1) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * r, s, 1)
    }


    fn dy_dr(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), r, 1) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), r, 1) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * s, r, 1)
    }


    fn dy_ds(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), s, 1) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), s, 1) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), s, 1) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * r, s, 1)
    }


    // fn jacobian(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
    //     r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
    //     tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    pub fn jacobian(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let unique_elements_of_rotation_matrix = 
            extract_unique_elements_of_rotation_matrix::<T, V>(ref_rotation_matrix)?;
        let shrinked_rotation_matrix = 
            ExtendedMatrix::create(T::from(3u8), T::from(3u8), 
            unique_elements_of_rotation_matrix, tolerance)?;

        let node_1 = ref_nodes.get(&node_1_number).unwrap();
        let node_2 = ref_nodes.get(&node_2_number).unwrap();
        let node_3 = ref_nodes.get(&node_3_number).unwrap();
        let node_4 = ref_nodes.get(&node_4_number).unwrap();
        let (node_1_x, node_1_y, node_1_z) = node_1.copy_coordinates();
        let (node_2_x, node_2_y, node_2_z) = node_2.copy_coordinates();
        let (node_3_x, node_3_y, node_3_z) = node_3.copy_coordinates(); 
        let (node_4_x, node_4_y, node_4_z) = node_4.copy_coordinates();

        let node_1_direction_vector = vec![
            node_1_x - node_3_x, node_1_y - node_3_y, node_1_z - node_3_z,
        ];
        let node_2_direction_vector = vec![
            node_2_x - node_3_x, node_2_y - node_3_y, node_2_z - node_3_z,
        ];
        let node_4_direction_vector = vec![
            node_4_x - node_3_x, node_4_y - node_3_y, node_4_z - node_3_z,
        ];
        
        let node_1_direction = ExtendedMatrix::create(T::from(3u8),
            T::from(1u8), node_1_direction_vector, tolerance)?;
        let node_2_direction = ExtendedMatrix::create(T::from(3u8),
            T::from(1u8), node_2_direction_vector, tolerance)?;
        let node_4_direction = ExtendedMatrix::create(T::from(3u8),
            T::from(1u8), node_4_direction_vector, tolerance)?;

        let transformed_node_1_direction = shrinked_rotation_matrix.multiply_by_matrix(&node_1_direction)?;
        let transformed_node_2_direction = shrinked_rotation_matrix.multiply_by_matrix(&node_2_direction)?;
        let transformed_node_4_direction = shrinked_rotation_matrix.multiply_by_matrix(&node_4_direction)?;

        let transformed_node_1_direction_x =
            matrix_element_value_extractor(T::from(0u8), T::from(0u8), &transformed_node_1_direction)?;
        let transformed_node_1_direction_y =
            matrix_element_value_extractor(T::from(1u8), T::from(0u8), &transformed_node_1_direction)?;
        let transformed_node_1_direction_z =
            matrix_element_value_extractor(T::from(2u8), T::from(0u8), &transformed_node_1_direction)?;
        let transformed_node_2_direction_x =
            matrix_element_value_extractor(T::from(0u8), T::from(0u8), &transformed_node_2_direction)?;
        let transformed_node_2_direction_y =
            matrix_element_value_extractor(T::from(1u8), T::from(0u8), &transformed_node_2_direction)?;
        let transformed_node_2_direction_z =
            matrix_element_value_extractor(T::from(2u8), T::from(0u8), &transformed_node_2_direction)?;
        let transformed_node_4_direction_x =
            matrix_element_value_extractor(T::from(0u8), T::from(0u8), &transformed_node_4_direction)?;
        let transformed_node_4_direction_y =
            matrix_element_value_extractor(T::from(1u8), T::from(0u8), &transformed_node_4_direction)?;
        let transformed_node_4_direction_z =
            matrix_element_value_extractor(T::from(2u8), T::from(0u8), &transformed_node_4_direction)?;

        // if compare_with_tolerance(transformed_node_1_direction_z - 
        //     transformed_node_2_direction_z, tolerance) != V::from(0f32) || 
        //     compare_with_tolerance(transformed_node_1_direction_z - 
        //         transformed_node_4_direction_z, tolerance) != V::from(0f32) || 
        //     compare_with_tolerance(transformed_node_2_direction_z - 
        //         transformed_node_4_direction_z, tolerance) != V::from(0f32)
        // {
        //     return Err("Quad membrane element Jacobian calculation: Incorrect nodes directions transformation!".into());
        // }
        
        let jacobian_elements = vec![
            QuadFullMemAuxFunctions::<T, V>::dx_dr(
                transformed_node_1_direction_x, transformed_node_2_direction_x, 
                V::from(0f32), transformed_node_4_direction_x, r, s),
            QuadFullMemAuxFunctions::<T, V>::dy_dr(
                transformed_node_1_direction_y, transformed_node_2_direction_y, 
                V::from(0f32), transformed_node_4_direction_y, r, s),
            QuadFullMemAuxFunctions::<T, V>::dx_ds(
                transformed_node_1_direction_x, transformed_node_2_direction_x, 
                V::from(0f32), transformed_node_4_direction_x, r, s),
            QuadFullMemAuxFunctions::<T, V>::dy_ds(
                transformed_node_1_direction_y, transformed_node_2_direction_y, 
                V::from(0f32), transformed_node_4_direction_y, r, s),
        ];

        let jacobian = ExtendedMatrix::create(
            T::from(2u8), T::from(2u8), jacobian_elements, tolerance)?;

        Ok(jacobian)
    }


    fn inverse_jacobian(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let jacobian = QuadFullMemAuxFunctions::<T, V>::jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)?;
        let inverse_jacobian = jacobian.inverse()?;
        Ok(inverse_jacobian)
    }


    fn determinant_of_jacobian(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<V, String>
    {
        let jacobian = QuadFullMemAuxFunctions::<T, V>::jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)?;
        let determinant_of_jacobian = jacobian.determinant_2x2()?;
        Ok(determinant_of_jacobian)
    }


    fn dh1_dr(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) + 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), r, 1) + 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh2_dr(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) - 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), r, 1) - 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh3_dr(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) - 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), r, 1) + 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh4_dr(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) + 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), r, 1) - 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh1_ds(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), s, 1) + 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) + 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh2_ds(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), s, 1) - 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) - 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh3_ds(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), s, 1) - 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) + 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh4_ds(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), s, 1) + 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) - 
        QuadFullMemAuxFunctions::<T, V>::derivative_x(
            QuadFullMemAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh_dx_dh_dy(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let inverse_jacobian = QuadFullMemAuxFunctions::inverse_jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, 
            ref_rotation_matrix, tolerance)?;
        let dh_dr_dh_ds = ExtendedMatrix::create(
            T::from(2u8), T::from(4u8), 
            vec![
                QuadFullMemAuxFunctions::<T, V>::dh1_dr(r, s), QuadFullMemAuxFunctions::<T, V>::dh2_dr(r, s), 
                QuadFullMemAuxFunctions::<T, V>::dh3_dr(r, s), QuadFullMemAuxFunctions::<T, V>::dh4_dr(r, s),
                QuadFullMemAuxFunctions::<T, V>::dh1_ds(r, s), QuadFullMemAuxFunctions::<T, V>::dh2_ds(r, s), 
                QuadFullMemAuxFunctions::<T, V>::dh3_ds(r, s), QuadFullMemAuxFunctions::<T, V>::dh4_ds(r, s),
            ], 
            tolerance)?;
        Ok(inverse_jacobian.multiply_by_matrix(&dh_dr_dh_ds)?)

    }


    pub fn strain_displacement_matrix(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let dh_dx_dh_dy_matrix = QuadFullMemAuxFunctions::<T, V>::dh_dx_dh_dy(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)?;

        let dh1_dx = matrix_element_value_extractor(T::from(0u8), T::from(0u8), &dh_dx_dh_dy_matrix)?;
        let dh2_dx = matrix_element_value_extractor(T::from(0u8), T::from(1u8), &dh_dx_dh_dy_matrix)?;
        let dh3_dx = matrix_element_value_extractor(T::from(0u8), T::from(2u8), &dh_dx_dh_dy_matrix)?;
        let dh4_dx = matrix_element_value_extractor(T::from(0u8), T::from(3u8), &dh_dx_dh_dy_matrix)?;

        let dh1_dy = matrix_element_value_extractor(T::from(1u8), T::from(0u8), &dh_dx_dh_dy_matrix)?;
        let dh2_dy = matrix_element_value_extractor(T::from(1u8), T::from(1u8), &dh_dx_dh_dy_matrix)?;
        let dh3_dy = matrix_element_value_extractor(T::from(1u8), T::from(2u8), &dh_dx_dh_dy_matrix)?;
        let dh4_dy = matrix_element_value_extractor(T::from(1u8), T::from(3u8), &dh_dx_dh_dy_matrix)?;

        let elements = vec![
            dh1_dx, V::from(0f32), V::from(0f32), dh2_dx, V::from(0f32), V::from(0f32), 
            dh3_dx, V::from(0f32), V::from(0f32), dh4_dx, V::from(0f32), V::from(0f32), 

            V::from(0f32), dh1_dy, V::from(0f32), V::from(0f32), dh2_dy, V::from(0f32),
            V::from(0f32), dh3_dy, V::from(0f32), V::from(0f32), dh4_dy, V::from(0f32),

            dh1_dy, dh1_dx, V::from(0f32), dh2_dy, dh2_dx, V::from(0f32), 
            dh3_dy, dh3_dx, V::from(0f32), dh4_dy, dh4_dx, V::from(0f32),
        ];

        let matrix = ExtendedMatrix::create(
            T::from(3u8), 
            QuadFullMemAuxFunctions::<T, V>::nodes_number() * QuadFullMemAuxFunctions::<T, V>::node_dof(), 
            elements, tolerance)?;

        Ok(matrix)
    }


    pub fn local_stiffness_matrix(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        young_modulus: V, poisson_ratio: V, thickness: V, alpha: V, r: V, s: V, 
        ref_local_stiffness_matrix: &ExtendedMatrix<T, V>, ref_nodes: &HashMap<T, FENode<V>>, 
        ref_rotation_matrix: &ExtendedMatrix<T, V>, tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let c_matrix_multiplier = young_modulus * thickness / (V::from(1f32) - poisson_ratio.my_powi(2));
        let mut c_matrix = ExtendedMatrix::create(
            T::from(3u8), T::from(3u8), 
            vec![
                V::from(1f32), poisson_ratio, V::from(0f32),
                poisson_ratio, V::from(1f32), V::from(0f32),
                V::from(0f32), V::from(0f32), (V::from(1f32) - poisson_ratio) / V::from(2f32),
            ], 
            tolerance)?;
        c_matrix.multiply_by_number(c_matrix_multiplier);

        let mut lhs_matrix = QuadFullMemAuxFunctions::strain_displacement_matrix(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, 
            ref_rotation_matrix, tolerance)?;
        lhs_matrix.transpose();

        let rhs_matrix = QuadFullMemAuxFunctions::strain_displacement_matrix(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, 
            ref_rotation_matrix, tolerance)?;

        return match (lhs_matrix.multiply_by_matrix(&c_matrix)?).multiply_by_matrix(&rhs_matrix)
        {
            Ok(mut matrix) =>
                {
                    matrix.multiply_by_number(QuadFullMemAuxFunctions::determinant_of_jacobian(
                        node_1_number, node_2_number, node_3_number, node_4_number, r, s, 
                        ref_nodes, ref_rotation_matrix, tolerance)? * alpha);

                    match ref_local_stiffness_matrix.add_matrix(&matrix)
                    {
                        Ok(matrix) => Ok(matrix),
                        Err(e) =>
                            Err(format!("Quad mem aux functions: Local stiffness matrix cannot be \
                                calculated! Reason: {}", e)),
                    }
                },
            Err(e) => Err(format!("Quad mem aux functions: Local stiffness matrix cannot be \
                                calculated! Reason: {}", e)),
        }
    }


    pub fn compose_node_dof_parameters<'a>(node_number: T)
        -> Result<Vec<DOFParameterData<T>>, &'a str>
    {
        let mut node_dof_parameters = Vec::new();
        for dof in 0..MEMBRANE_NODE_DOF
        {
            let dof_parameter =
                GlobalDOFParameter::iterator().nth(dof)
                    .ok_or("Quad mem aux functions: Could not find dof parameter!")?;
            let dof_parameter = DOFParameterData::create(
                node_number, *dof_parameter);
            node_dof_parameters.push(dof_parameter);
        }
        Ok(node_dof_parameters)
    }


    pub fn extract_column_matrix_values(column_matrix: &ExtendedMatrix<T, V>)
        -> Result<Vec<V>, String>
    {
        let mut values = Vec::new();
        let shape = column_matrix.copy_shape();

        let mut row = T::from(0u8);
        while row < shape.0
        {
            let mut column = T::from(0u8);
            while column < shape.1
            {
                let value = matrix_element_value_extractor(row, column, column_matrix)?;
                values.push(value);
                column += T::from(1u8);
            }
            row += T::from(1u8);
        }
        Ok(values)
    }
}
