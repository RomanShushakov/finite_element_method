use std::ops::{Add, Sub, Div, Rem, SubAssign, Mul, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::collections::HashMap;

use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::matrix_element_position::MatrixElementPosition;

use crate::my_float::MyFloatTrait;

use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData, GlobalDOFParameter};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::membrane::consts::{MEM4N4IP_NODES_NUMBER, MEMBRANE_NODE_DOF};
use crate::fem::finite_elements::functions::
{
    compare_with_tolerance, extract_unique_elements_of_rotation_matrix
};


pub struct QuadMemAuxFunctions<T, V>(T, V);


impl<T, V> QuadMemAuxFunctions<T, V>
    where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
             From<u8> + Ord + 'static,
          V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + From<f32> +
             Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
             MyFloatTrait + PartialOrd + 'static,
{
    fn length(node_1_number: T, node_2_number: T, nodes: &HashMap<T, FENode<V>>) -> V
    {
        let node_1 = nodes.get(&node_1_number).unwrap();
        let node_2 = nodes.get(&node_2_number).unwrap();

        ((node_1.copy_x() - node_2.copy_x()).my_powi(2) + (node_1.copy_y() - node_2.copy_y()).my_powi(2) +
        (node_1.copy_z() - node_2.copy_z()).my_powi(2)).my_sqrt()
    }


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
        tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let node_2 = nodes.get(&node_2_number).unwrap();
        let node_3 = nodes.get(&node_3_number).unwrap();
        let node_4 = nodes.get(&node_4_number).unwrap();

        let edge_3_4_x = node_4.copy_x() - node_3.copy_x();
        let edge_3_4_y = node_4.copy_y() - node_3.copy_y();
        let edge_3_4_z = node_4.copy_z() - node_3.copy_z();

        let edge_3_4_length = QuadMemAuxFunctions::<T, V>::length(node_3_number, node_4_number, nodes);

        let (u, v, w) = (edge_3_4_length, V::from(0f32), V::from(0f32));
        let alpha = ((edge_3_4_x * u + edge_3_4_y * v + edge_3_4_z * w) / 
            (edge_3_4_length * edge_3_4_length)).my_acos();
        let (rotation_axis_coord_x, mut rotation_axis_coord_y,
            mut rotation_axis_coord_z) = (V::from(0f32), V::from(0f32), V::from(0f32));
        if edge_3_4_x != V::from(0f32) && edge_3_4_y == V::from(0f32) && edge_3_4_z == V::from(0f32)
        {
            rotation_axis_coord_z = edge_3_4_x;
        }
        else
        {
            rotation_axis_coord_y = edge_3_4_z * edge_3_4_length;
            rotation_axis_coord_z = V::from(-1f32) * edge_3_4_y * edge_3_4_length;
        }
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

        let interim_rotation_matrix = ExtendedMatrix::create(
            3, 3,
            vec![q_11, q_12, q_13, q_21, q_22, q_23, q_31, q_32, q_33],
            tolerance)?;

        let edge_3_2_x = node_2.copy_x() - node_3.copy_x();
        let edge_3_2_y = node_2.copy_y() - node_3.copy_y();
        let edge_3_2_z = node_2.copy_z() - node_3.copy_z();

        let normal_through_node_3_x = edge_3_4_y * edge_3_2_z - edge_3_4_z * edge_3_2_y;
        let normal_through_node_3_y = edge_3_4_z * edge_3_2_x - edge_3_4_x * edge_3_2_z;
        let normal_through_node_3_z = edge_3_4_x * edge_3_2_y - edge_3_4_y * edge_3_2_x;

        let normal_through_node_3 = ExtendedMatrix::create(
            3, 1,
            vec![normal_through_node_3_x, 
                normal_through_node_3_y, normal_through_node_3_z],
            tolerance)?;

        let transformed_normal_through_node_3 =
            interim_rotation_matrix.multiply_by_matrix(&normal_through_node_3)?;

        let transformed_normal_through_node_3_x =
            transformed_normal_through_node_3.copy_element_value_or_zero(
                MatrixElementPosition::create(0, 0))?;

        let transformed_normal_through_node_3_y =
            transformed_normal_through_node_3.copy_element_value_or_zero(
                MatrixElementPosition::create(1, 0))?;

        let transformed_normal_through_node_3_z =
            transformed_normal_through_node_3.copy_element_value_or_zero(
                MatrixElementPosition::create(2, 0))?;

        let angle_between_normal_direction_and_axis_t =
            (transformed_normal_through_node_3_z /
            (transformed_normal_through_node_3_x.my_powi(2) +
            transformed_normal_through_node_3_y.my_powi(2) +
            transformed_normal_through_node_3_z.my_powi(2))
                .my_sqrt()
            ).my_acos();
        
        let c_x = compare_with_tolerance(edge_3_4_x / edge_3_4_length, tolerance);
        let c_y = compare_with_tolerance(edge_3_4_y / edge_3_4_length, tolerance);
        let c_z = compare_with_tolerance(edge_3_4_z / edge_3_4_length, tolerance);
        let c_xz = compare_with_tolerance(
            (c_x.my_powi(2) + c_z.my_powi(2)).my_sqrt(), tolerance);

        let c = angle_between_normal_direction_and_axis_t.my_cos();
        let s = angle_between_normal_direction_and_axis_t.my_sin();

        let r_11 = if c_xz != V::from(0f32) { c_x } else { V::from(0f32) };
        let r_12 = c_y;
        let r_13 = if c_xz != V::from(0f32) { c_z } else { V::from(0f32) };
        let r_21 = if c_xz != V::from(0f32) { (V::from(-1f32) * c_x * c_y * c - c_z * s) / c_xz }
            else { V::from(-1f32) * c_y * c };
        let r_22 = if c_xz != V::from(0f32) { c_xz * c } else { V::from(0f32) };
        let r_23 = if c_xz != V::from(0f32)
            { (V::from(-1f32) * c_y * c_z * c + c_x * s) / c_xz } else { s };
        let r_31 = if c_xz != V::from(0f32)
            { (c_x * c_y * s - c_z * c) / c_xz } else { c_y * s };
        let r_32 = if c_xz != V::from(0f32)
            { V::from(-1f32) * c_xz * s } else { V::from(0f32) };
        let r_33 = if c_xz != V::from(0f32)
            { (c_y * c_z * s + c_x * c) / c_xz } else { c };

        let rotation_matrix = ExtendedMatrix::create(
            QuadMemAuxFunctions::<T, V>::nodes_number() * QuadMemAuxFunctions::<T, V>::node_dof(),
            QuadMemAuxFunctions::<T, V>::nodes_number() * QuadMemAuxFunctions::<T, V>::node_dof(),
            vec![
                [r_11, r_12, r_13], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [r_21, r_22, r_23], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [r_31, r_32, r_33], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],

                [V::from(0f32); MEMBRANE_NODE_DOF], [r_11, r_12, r_13],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [r_21, r_22, r_23],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [r_31, r_32, r_33],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],

                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [r_11, r_12, r_13], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [r_21, r_22, r_23], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [r_31, r_32, r_33], [V::from(0f32); MEMBRANE_NODE_DOF],

                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [r_11, r_12, r_13],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [r_21, r_22, r_23],
                [V::from(0f32); MEMBRANE_NODE_DOF], [V::from(0f32); MEMBRANE_NODE_DOF],
                [V::from(0f32); MEMBRANE_NODE_DOF], [r_31, r_32, r_33],
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
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * s, r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), r, 1) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * s, r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * s, r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), r, 1) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * s, r, 1)
    }


    fn dx_ds(x_1: V, x_2: V, x_3: V, x_4: V, r: V, s: V) -> V
    {
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * r, s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), s, 1) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * r, s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), s, 1) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * r, s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * r, s, 1)
    }


    fn dy_dr(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
    {
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * s, r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), r, 1) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * s, r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * s, r, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), r, 1) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * s, r, 1)
    }


    fn dy_ds(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
    {
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * r, s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), s, 1) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * r, s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), s, 1) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * r, s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), s, 1) +
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadMemAuxFunctions::<T, V>::derivative_x(
            QuadMemAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * r, s, 1)
    }


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

        let transformed_transformed_node_1_direction_x =
            transformed_node_1_direction.copy_element_value_or_zero(
                MatrixElementPosition::create(T::from(0u8), T::from(0u8)))?;
        let transformed_transformed_node_1_direction_y =
            transformed_node_1_direction.copy_element_value_or_zero(
                MatrixElementPosition::create(T::from(1u8), T::from(0u8)))?;
        let transformed_transformed_node_1_direction_z =
            transformed_node_1_direction.copy_element_value_or_zero(
                MatrixElementPosition::create(T::from(2u8), T::from(0u8)))?;
        let transformed_transformed_node_2_direction_x =
            transformed_node_2_direction.copy_element_value_or_zero(
                MatrixElementPosition::create(T::from(0u8), T::from(0u8)))?;
        let transformed_transformed_node_2_direction_y =
            transformed_node_2_direction.copy_element_value_or_zero(
                MatrixElementPosition::create(T::from(1u8), T::from(0u8)))?;
        let transformed_transformed_node_2_direction_z =
            transformed_node_2_direction.copy_element_value_or_zero(
                MatrixElementPosition::create(T::from(2u8), T::from(0u8)))?;
        let transformed_transformed_node_4_direction_x =
            transformed_node_4_direction.copy_element_value_or_zero(
                MatrixElementPosition::create(T::from(0u8), T::from(0u8)))?;
        let transformed_transformed_node_4_direction_y =
            transformed_node_4_direction.copy_element_value_or_zero(
                MatrixElementPosition::create(T::from(1u8), T::from(0u8)))?;
        let transformed_transformed_node_4_direction_z =
            transformed_node_4_direction.copy_element_value_or_zero(
                MatrixElementPosition::create(T::from(2u8), T::from(0u8)))?;

        println!("{:?}, {:?}, {:?}", transformed_transformed_node_1_direction_z,
            transformed_transformed_node_2_direction_z, transformed_transformed_node_4_direction_z);
        
        let jacobian_elements = vec![
            QuadMemAuxFunctions::<T, V>::dx_dr(
                transformed_transformed_node_1_direction_x, transformed_transformed_node_2_direction_x, 
                V::from(0f32), transformed_transformed_node_4_direction_x, r, s),
            QuadMemAuxFunctions::<T, V>::dy_dr(
                transformed_transformed_node_1_direction_y, transformed_transformed_node_2_direction_y, 
                V::from(0f32), transformed_transformed_node_4_direction_y, r, s),
            QuadMemAuxFunctions::<T, V>::dx_ds(
                transformed_transformed_node_1_direction_x, transformed_transformed_node_2_direction_x, 
                V::from(0f32), transformed_transformed_node_4_direction_x, r, s),
            QuadMemAuxFunctions::<T, V>::dy_ds(
                transformed_transformed_node_1_direction_y, transformed_transformed_node_2_direction_y, 
                V::from(0f32), transformed_transformed_node_4_direction_y, r, s),
        ];

        let jacobian = ExtendedMatrix::create(
            T::from(2u8), T::from(2u8), jacobian_elements, tolerance)?;

        Ok(jacobian)
    }


    // fn inverse_jacobian(node_1_number: T, node_2_number: T, r: V, nodes: &HashMap<T, FENode<V>>)
    //     -> V
    // {
    //     V::from(1f32) / TrussAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    // }


    // fn determinant_of_jacobian(node_1_number: T, node_2_number: T, r: V,
    //     nodes: &HashMap<T, FENode<V>>) -> V
    // {
    //     TrussAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    // }


    // fn dh1_dr(r: V) -> V
    // {
    //     TrussAuxFunctions::<T, V>::derivative_x(
    //         TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) -
    //     TrussAuxFunctions::<T, V>::derivative_x(
    //         TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), r, 1)
    // }


    // fn dh2_dr(r: V) -> V
    // {
    //     TrussAuxFunctions::<T, V>::derivative_x(
    //         TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) +
    //     TrussAuxFunctions::<T, V>::derivative_x(
    //         TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), r, 1)
    // }


    // pub fn strain_displacement_matrix(node_1_number: T, node_2_number: T, r: V, tolerance: V,
    //     nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    // {
    //     let elements = vec![TrussAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32),
    //         V::from(0f32), TrussAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32), V::from(0f32)];
    //     let mut matrix = ExtendedMatrix::create(T::from(1u8),
    //         TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
    //         elements, tolerance)?;
    //     let inverse_jacobian = TrussAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
    //         r, nodes);
    //     matrix.multiply_by_number(inverse_jacobian);
    //     Ok(matrix)
    // }


    // pub fn area(area_1: V, area_2: Option<V>, r: V) -> V
    // {
    //     if let Some(area_2) = area_2
    //     {
    //         (area_2 - area_1) / V::from(2f32) * r + area_1 -
    //             (area_2 - area_1) / V::from(2f32) * V::from(-1f32)
    //     }
    //     else
    //     {
    //         area_1
    //     }
    // }


    // pub fn local_stiffness_matrix(node_1_number: T, node_2_number: T, young_modulus: V, area_1: V,
    //     area_2: Option<V>, alpha: V, r: V, local_stiffness_matrix: &ExtendedMatrix<T, V>,
    //     tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    // {
    //     let current_area = TrussAuxFunctions::<T, V>::area(area_1, area_2, r);

    //     let mut lhs_matrix = TrussAuxFunctions::strain_displacement_matrix(
    //         node_1_number, node_2_number, r, tolerance, nodes)?;

    //     lhs_matrix.transpose();

    //     lhs_matrix.multiply_by_number(young_modulus * current_area);

    //     let rhs_matrix = TrussAuxFunctions::strain_displacement_matrix(
    //             node_1_number, node_2_number, r, tolerance, nodes)?;

    //     return match lhs_matrix.multiply_by_matrix(&rhs_matrix)
    //     {
    //         Ok(mut matrix) =>
    //             {
    //                 matrix.multiply_by_number(TrussAuxFunctions::determinant_of_jacobian(
    //                     node_1_number, node_2_number, r, nodes) * alpha);

    //                 match local_stiffness_matrix.add_matrix(&matrix)
    //                 {
    //                     Ok(matrix) => Ok(matrix),
    //                     Err(e) =>
    //                         Err(format!("Truss2n2ip: Local stiffness matrix cannot be \
    //                             calculated! Reason: {}", e)),
    //                 }
    //             },
    //         Err(e) => Err(format!("Truss2n2ip: Local stiffness matrix cannot be \
    //                             calculated! Reason: {}", e)),
    //     }
    // }


    // pub fn compose_node_dof_parameters<'a>(node_number: T)
    //     -> Result<Vec<DOFParameterData<T>>, &'a str>
    // {
    //     let mut node_dof_parameters = Vec::new();
    //     for dof in 0..TRUSS_NODE_DOF
    //     {
    //         let dof_parameter =
    //             GlobalDOFParameter::iterator().nth(dof)
    //                 .ok_or("Truss2n2ip: Could not find dof parameter!")?;
    //         let dof_parameter = DOFParameterData::create(
    //             node_number, *dof_parameter);
    //         node_dof_parameters.push(dof_parameter);
    //     }
    //     Ok(node_dof_parameters)
    // }


    // pub fn extract_column_matrix_values(column_matrix: &ExtendedMatrix<T, V>)
    //     -> Result<Vec<V>, String>
    // {
    //     let mut values = Vec::new();
    //     let shape = column_matrix.copy_shape();

    //     let mut row = T::from(0u8);
    //     while row < shape.0
    //     {
    //         let mut column = T::from(0u8);
    //         while column < shape.1
    //         {
    //             let value = column_matrix.copy_element_value_or_zero(
    //                 MatrixElementPosition::create(row, column))?;;
    //             values.push(value);
    //             column += T::from(1u8);
    //         }
    //         row += T::from(1u8);
    //     }
    //     Ok(values)
    // }
}
