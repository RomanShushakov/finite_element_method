use std::ops::{Add, Sub, Div, Rem, SubAssign, Mul, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::collections::HashMap;

use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::matrix_element_value_extractor;

use extended_matrix_float::MyFloatTrait;

use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData, GlobalDOFParameter};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::plate::consts::{PLATE4N4IP_NODES_NUMBER, PLATE_NODE_DOF};
use crate::fem::finite_elements::functions::
{
    compare_with_tolerance, extract_unique_elements_of_rotation_matrix
};


pub struct QuadFullPlateAuxFunctions<T, V>(T, V);


impl<T, V> QuadFullPlateAuxFunctions<T, V>
    where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
             From<u8> + Ord + 'static,
          V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + From<f32> +
             Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
             MyFloatTrait + PartialOrd + 'static,
{
    pub fn nodes_number() -> T
    {
        let mut n = T::from(0u8);
        (0..PLATE4N4IP_NODES_NUMBER).for_each(|_| n += T::from(1u8));
        n
    }


    pub fn node_dof() -> T
    {
        let mut n = T::from(0u8);
        (0..PLATE_NODE_DOF).for_each(|_| n += T::from(1u8));
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
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(),
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(),
            vec![
                [q_11, q_12, q_13, V::from(0f32), V::from(0f32), V::from(0f32)],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF],

                [q_21, q_22, q_23, V::from(0f32), V::from(0f32), V::from(0f32)], 
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],

                [q_31, q_32, q_33, V::from(0f32), V::from(0f32), V::from(0f32)], 
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32), V::from(0f32), V::from(0f32), q_11, q_12, q_13],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32), V::from(0f32), V::from(0f32), q_21, q_22, q_23],
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32), V::from(0f32), V::from(0f32), q_31, q_32, q_33],
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],


                [V::from(0f32); PLATE_NODE_DOF], 
                [q_11, q_12, q_13, V::from(0f32), V::from(0f32), V::from(0f32)],
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF],
                [q_21, q_22, q_23, V::from(0f32), V::from(0f32), V::from(0f32)], 
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF],
                [q_31, q_32, q_33, V::from(0f32), V::from(0f32), V::from(0f32)], 
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32), V::from(0f32), V::from(0f32), q_11, q_12, q_13],
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32), V::from(0f32), V::from(0f32), q_21, q_22, q_23],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32), V::from(0f32), V::from(0f32), q_31, q_32, q_33],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],


                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [q_11, q_12, q_13, V::from(0f32), V::from(0f32), V::from(0f32)],
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [q_21, q_22, q_23, V::from(0f32), V::from(0f32), V::from(0f32)], 
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [q_31, q_32, q_33, V::from(0f32), V::from(0f32), V::from(0f32)], 
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32), V::from(0f32), V::from(0f32), q_11, q_12, q_13],
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32), V::from(0f32), V::from(0f32), q_21, q_22, q_23],
                [V::from(0f32); PLATE_NODE_DOF],

                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32), V::from(0f32), V::from(0f32), q_31, q_32, q_33],
                [V::from(0f32); PLATE_NODE_DOF],


                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [q_11, q_12, q_13, V::from(0f32), V::from(0f32), V::from(0f32)],

                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [q_21, q_22, q_23, V::from(0f32), V::from(0f32), V::from(0f32)], 

                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [q_31, q_32, q_33, V::from(0f32), V::from(0f32), V::from(0f32)], 

                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32), V::from(0f32), V::from(0f32), q_11, q_12, q_13],
                
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32), V::from(0f32), V::from(0f32), q_21, q_22, q_23],
                
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32); PLATE_NODE_DOF], 
                [V::from(0f32); PLATE_NODE_DOF],
                [V::from(0f32), V::from(0f32), V::from(0f32), q_31, q_32, q_33],
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
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * s, r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), r, 1) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * s, r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * s, r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), r, 1) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * s, r, 1)
    }


    fn dx_ds(x_1: V, x_2: V, x_3: V, x_4: V, r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32), s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.25f32) * r, s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32), s, 1) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.25f32) * r, s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32), s, 1) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_3 * V::from(0.25f32) * r, s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32), s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, x_4 * V::from(0.25f32) * r, s, 1)
    }


    fn dy_dr(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * s, r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), r, 1) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * s, r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * s, r, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), r, 1) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * s, r, 1)
    }


    fn dy_ds(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32), s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_1 * V::from(0.25f32) * r, s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32), s, 1) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_2 * V::from(0.25f32) * r, s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32), s, 1) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_3 * V::from(0.25f32) * r, s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32), s, 1) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, y_4 * V::from(0.25f32) * r, s, 1)
    }


    fn extract_transformed_directions_of_nodes(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>, tolerance: V) 
        -> Result<Vec<(V, V, V)>, String>
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

        if transformed_node_1_direction_z != transformed_node_2_direction_z || 
            transformed_node_1_direction_z != transformed_node_4_direction_z ||
            transformed_node_2_direction_z != transformed_node_4_direction_z
        {
            return Err("Quad plate element: Transformed directions of nodes calculation: /
                Incorrect nodes directions transformation!".into());
        }

        let transformed_directions_of_nodes = vec![
            (transformed_node_1_direction_x, transformed_node_1_direction_y, transformed_node_1_direction_z),
            (transformed_node_2_direction_x, transformed_node_2_direction_y, transformed_node_2_direction_z),
            (transformed_node_4_direction_x, transformed_node_4_direction_y, transformed_node_4_direction_z),
        ];
        Ok(transformed_directions_of_nodes)
    }


    pub fn jacobian(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let transformed_directions_of_nodes = 
            QuadFullPlateAuxFunctions::extract_transformed_directions_of_nodes(
                node_1_number, node_2_number, node_3_number, node_4_number, ref_nodes, ref_rotation_matrix, tolerance)?;

        if transformed_directions_of_nodes.len() != 3
        {
            return Err("Quad plate element Jacobian calculation: Incorrect quantity of transformed directions of nodes!".into());
        }

        let transformed_node_1_direction_x = transformed_directions_of_nodes[0].0;
        let transformed_node_1_direction_y = transformed_directions_of_nodes[0].1;
        let transformed_node_2_direction_x = transformed_directions_of_nodes[1].0;
        let transformed_node_2_direction_y = transformed_directions_of_nodes[1].1;
        let transformed_node_4_direction_x = transformed_directions_of_nodes[2].0;
        let transformed_node_4_direction_y = transformed_directions_of_nodes[2].1;
        
        let jacobian_elements = vec![
            QuadFullPlateAuxFunctions::<T, V>::dx_dr(
                transformed_node_1_direction_x, transformed_node_2_direction_x, 
                V::from(0f32), transformed_node_4_direction_x, r, s),
            QuadFullPlateAuxFunctions::<T, V>::dy_dr(
                transformed_node_1_direction_y, transformed_node_2_direction_y, 
                V::from(0f32), transformed_node_4_direction_y, r, s),
            QuadFullPlateAuxFunctions::<T, V>::dx_ds(
                transformed_node_1_direction_x, transformed_node_2_direction_x, 
                V::from(0f32), transformed_node_4_direction_x, r, s),
            QuadFullPlateAuxFunctions::<T, V>::dy_ds(
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
        let jacobian = QuadFullPlateAuxFunctions::<T, V>::jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)?;
        let inverse_jacobian = jacobian.inverse()?;
        Ok(inverse_jacobian)
    }


    fn determinant_of_jacobian(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<V, String>
    {
        let jacobian = QuadFullPlateAuxFunctions::<T, V>::jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)?;
        let determinant_of_jacobian = jacobian.determinant_2x2()?;
        Ok(determinant_of_jacobian)
    }


    fn dh1_dr(r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) + 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), r, 1) + 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh2_dr(r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) - 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), r, 1) - 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh3_dr(r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) - 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), r, 1) + 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh4_dr(r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) + 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), r, 1) - 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh1_ds(r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), s, 1) + 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) + 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh2_ds(r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), s, 1) - 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) - 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh3_ds(r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), s, 1) - 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) + 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh4_ds(r: V, s: V) -> V
    {
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32), s, 1) + 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) - 
        QuadFullPlateAuxFunctions::<T, V>::derivative_x(
            QuadFullPlateAuxFunctions::<T, V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh_dx_dh_dy(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let inverse_jacobian = QuadFullPlateAuxFunctions::inverse_jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, 
            ref_rotation_matrix, tolerance)?;
        let dh_dr_dh_ds = ExtendedMatrix::create(
            T::from(2u8), T::from(4u8), 
            vec![
                QuadFullPlateAuxFunctions::<T, V>::dh1_dr(r, s), QuadFullPlateAuxFunctions::<T, V>::dh2_dr(r, s), 
                QuadFullPlateAuxFunctions::<T, V>::dh3_dr(r, s), QuadFullPlateAuxFunctions::<T, V>::dh4_dr(r, s),
                QuadFullPlateAuxFunctions::<T, V>::dh1_ds(r, s), QuadFullPlateAuxFunctions::<T, V>::dh2_ds(r, s), 
                QuadFullPlateAuxFunctions::<T, V>::dh3_ds(r, s), QuadFullPlateAuxFunctions::<T, V>::dh4_ds(r, s),
            ], 
            tolerance)?;
        Ok(inverse_jacobian.multiply_by_matrix(&dh_dr_dh_ds)?)

    }


    fn h1_r_s(r: V, s: V) -> V
    {
        V::from(0.25f32) * (V::from(1f32) + r) * (V::from(1f32) + s)
    }


    fn h2_r_s(r: V, s: V) -> V
    {
        V::from(0.25f32) * (V::from(1f32) - r) * (V::from(1f32) + s)
    }


    fn h3_r_s(r: V, s: V) -> V
    {
        V::from(0.25f32) * (V::from(1f32) - r) * (V::from(1f32) - s)
    }


    fn h4_r_s(r: V, s: V) -> V
    {
        V::from(0.25f32) * (V::from(1f32) + r) * (V::from(1f32) - s)
    }


    pub fn strain_displacement_matrix_mem(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let dh_dx_dh_dy_matrix = QuadFullPlateAuxFunctions::<T, V>::dh_dx_dh_dy(
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
            dh1_dx, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            dh2_dx, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            dh3_dx, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            dh4_dx, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),

            V::from(0f32), dh1_dy, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), dh2_dy, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), dh3_dy, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), dh4_dy, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),

            dh1_dy, dh1_dx, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            dh2_dy, dh2_dx, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            dh3_dy, dh3_dx, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            dh4_dy, dh4_dx, V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
        ];

        let matrix = ExtendedMatrix::create(
            T::from(3u8), 
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(), 
            elements, tolerance)?;

        Ok(matrix)
    }


    pub fn strain_displacement_matrix_plate_bending(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let dh_dx_dh_dy_matrix = QuadFullPlateAuxFunctions::<T, V>::dh_dx_dh_dy(
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
            V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(-1f32) * dh1_dx, V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(-1f32) * dh2_dx, V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(-1f32) * dh3_dx, V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(-1f32) * dh4_dx, V::from(0f32),

            V::from(0f32), V::from(0f32), V::from(0f32), dh1_dy, V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), dh2_dy, V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), dh3_dy, V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), dh4_dy, V::from(0f32), V::from(0f32),

            V::from(0f32), V::from(0f32), V::from(0f32), dh1_dx, V::from(-1f32) * dh1_dy, V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), dh2_dx, V::from(-1f32) * dh2_dy, V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), dh3_dx, V::from(-1f32) * dh3_dy, V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), dh4_dx, V::from(-1f32) * dh4_dy, V::from(0f32),
        ];

        let matrix = ExtendedMatrix::create(
            T::from(3u8), 
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(), 
            elements, tolerance)?;

        Ok(matrix)
    }


    // pub fn strain_displacement_matrix_plate_shear(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
    //     r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
    //     tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    // {
    //     let dh_dx_dh_dy_matrix = QuadFullPlateAuxFunctions::<T, V>::dh_dx_dh_dy(
    //         node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)?;

    //     let dh1_dx = matrix_element_value_extractor(T::from(0u8), T::from(0u8), &dh_dx_dh_dy_matrix)?;
    //     let dh2_dx = matrix_element_value_extractor(T::from(0u8), T::from(1u8), &dh_dx_dh_dy_matrix)?;
    //     let dh3_dx = matrix_element_value_extractor(T::from(0u8), T::from(2u8), &dh_dx_dh_dy_matrix)?;
    //     let dh4_dx = matrix_element_value_extractor(T::from(0u8), T::from(3u8), &dh_dx_dh_dy_matrix)?;

    //     let dh1_dy = matrix_element_value_extractor(T::from(1u8), T::from(0u8), &dh_dx_dh_dy_matrix)?;
    //     let dh2_dy = matrix_element_value_extractor(T::from(1u8), T::from(1u8), &dh_dx_dh_dy_matrix)?;
    //     let dh3_dy = matrix_element_value_extractor(T::from(1u8), T::from(2u8), &dh_dx_dh_dy_matrix)?;
    //     let dh4_dy = matrix_element_value_extractor(T::from(1u8), T::from(3u8), &dh_dx_dh_dy_matrix)?;

    //     let h_1 = QuadFullPlateAuxFunctions::<T, V>::h1_r_s(r, s);
    //     let h_2 = QuadFullPlateAuxFunctions::<T, V>::h2_r_s(r, s);
    //     let h_3 = QuadFullPlateAuxFunctions::<T, V>::h3_r_s(r, s);
    //     let h_4 = QuadFullPlateAuxFunctions::<T, V>::h4_r_s(r, s);

    //     let elements = vec![
    //         V::from(0f32), V::from(0f32), dh1_dx, V::from(0f32), h_1, V::from(0f32),
    //         V::from(0f32), V::from(0f32), dh2_dx, V::from(0f32), h_2, V::from(0f32),
    //         V::from(0f32), V::from(0f32), dh3_dx, V::from(0f32), h_3, V::from(0f32),
    //         V::from(0f32), V::from(0f32), dh4_dx, V::from(0f32), h_4, V::from(0f32),

    //         V::from(0f32), V::from(0f32), dh1_dy, V::from(-1f32) * h_1, V::from(0f32), V::from(0f32),
    //         V::from(0f32), V::from(0f32), dh2_dy, V::from(-1f32) * h_2, V::from(0f32), V::from(0f32),
    //         V::from(0f32), V::from(0f32), dh3_dy, V::from(-1f32) * h_3, V::from(0f32), V::from(0f32),
    //         V::from(0f32), V::from(0f32), dh4_dy, V::from(-1f32) * h_4, V::from(0f32), V::from(0f32),
    //     ];

    //     let matrix = ExtendedMatrix::create(
    //         T::from(2u8), 
    //         QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(), 
    //         elements, tolerance)?;

    //     Ok(matrix)
    // }


    pub fn strain_displacement_matrix_plate_shear(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
        tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        let transformed_directions_of_nodes = 
            QuadFullPlateAuxFunctions::extract_transformed_directions_of_nodes(
                node_1_number, node_2_number, node_3_number, node_4_number, ref_nodes, ref_rotation_matrix, tolerance)?;

        if transformed_directions_of_nodes.len() != 3
        {
            return Err("Quad plate element strain displacement matrix for shear calculation: / 
                Incorrect quantity of transformed directions of nodes!".into());
        }

        let x_1 = transformed_directions_of_nodes[0].0;
        let y_1 = transformed_directions_of_nodes[0].1;
        let x_2 = transformed_directions_of_nodes[1].0;
        let y_2 = transformed_directions_of_nodes[1].1;
        let x_3 = V::from(0f32);
        let y_3 = V::from(0f32);
        let x_4 = transformed_directions_of_nodes[2].0;
        let y_4 = transformed_directions_of_nodes[2].1;

        let a_x = x_1 - x_2 - x_3 + x_4;
        let b_x = x_1 - x_2 + x_3 - x_4;
        let c_x = x_1 + x_2 - x_3 - x_4;
        let a_y = y_1 - y_2 - y_3 + y_4;
        let b_y = y_1 - y_2 + y_3 - y_4;
        let c_y = y_1 + y_2 - y_3 - y_4;

        let determinant_of_jacobian = QuadFullPlateAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)?;

        let gamma_rz_multiplier = ((c_x + r * b_x).my_powi(2) + (c_y + r * b_y).my_powi(2)).my_sqrt() / 
            (V::from(8f32) * determinant_of_jacobian);
        
        let gamma_sz_multiplier = ((a_x + s * b_x).my_powi(2) + (a_y + s * b_y).my_powi(2)).my_sqrt() / 
            (V::from(8f32) * determinant_of_jacobian);

        let elements = vec![
            V::from(0f32), 
            V::from(0f32), 
            (V::from(1f32) + s) / V::from(2f32) * gamma_rz_multiplier,
            (V::from(1f32) + s) * V::from(-1f32) * (y_1 - y_2) / V::from(4f32) * gamma_rz_multiplier,
            (V::from(1f32) + s) * (x_1 - x_2) / V::from(4f32) * gamma_rz_multiplier,
            V::from(0f32),

            V::from(0f32), 
            V::from(0f32), 
            V::from(-1f32) * (V::from(1f32) + s) / V::from(2f32) * gamma_rz_multiplier,
            (V::from(1f32) + s) * V::from(-1f32) * (y_1 - y_2) / V::from(4f32) * gamma_rz_multiplier,
            (V::from(1f32) + s) * (x_1 - x_2) / V::from(4f32) * gamma_rz_multiplier,
            V::from(0f32),

            V::from(0f32), 
            V::from(0f32), 
            V::from(-1f32) * (V::from(1f32) - s) / V::from(2f32) * gamma_rz_multiplier,
            (V::from(1f32) - s) * V::from(-1f32) * (y_4 - y_3) / V::from(4f32) * gamma_rz_multiplier,
            (V::from(1f32) - s) * (x_4 - x_3) / V::from(4f32) * gamma_rz_multiplier,
            V::from(0f32),

            V::from(0f32), 
            V::from(0f32), 
            (V::from(1f32) - s) / V::from(2f32) * gamma_rz_multiplier,
            (V::from(1f32) - s) * V::from(-1f32) * (y_4 - y_3) / V::from(4f32) * gamma_rz_multiplier,
            (V::from(1f32) - s) * (x_4 - x_3) / V::from(4f32) * gamma_rz_multiplier,
            V::from(0f32),


            V::from(0f32), 
            V::from(0f32), 
            (V::from(1f32) + r) / V::from(2f32) * gamma_sz_multiplier,
            (V::from(1f32) + r) * V::from(-1f32) * (y_1 - y_4) / V::from(4f32) * gamma_sz_multiplier,
            (V::from(1f32) + r) * (x_1 - x_4) / V::from(4f32) * gamma_sz_multiplier,
            V::from(0f32),

            V::from(0f32), 
            V::from(0f32), 
            (V::from(1f32) - r) / V::from(2f32) * gamma_sz_multiplier,
            (V::from(1f32) - r) * V::from(-1f32) * (y_2 - y_3) / V::from(4f32) * gamma_sz_multiplier,
            (V::from(1f32) - r) * (x_2 - x_3) / V::from(4f32) * gamma_sz_multiplier,
            V::from(0f32),

            V::from(0f32), 
            V::from(0f32), 
            V::from(-1f32) * (V::from(1f32) - r) / V::from(2f32) * gamma_sz_multiplier,
            (V::from(1f32) - r) * V::from(-1f32) * (y_2 - y_3) / V::from(4f32) * gamma_sz_multiplier,
            (V::from(1f32) - r) * (x_2 - x_3) / V::from(4f32) * gamma_sz_multiplier,
            V::from(0f32),

            V::from(0f32), 
            V::from(0f32), 
            V::from(-1f32) * (V::from(1f32) + r) / V::from(2f32) * gamma_sz_multiplier,
            (V::from(1f32) + r) * V::from(-1f32) * (y_1 - y_4) / V::from(4f32) * gamma_sz_multiplier,
            (V::from(1f32) + r) * (x_1 - x_4) / V::from(4f32) * gamma_sz_multiplier,
            V::from(0f32),
        ];

        let matrix = ExtendedMatrix::create(
            T::from(2u8), 
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(), 
            elements, tolerance)?;

        Ok(matrix)
    }


    pub fn local_stiffness_matrix(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
        young_modulus: V, poisson_ratio: V, thickness: V, shear_factor: V, alpha: V, r: V, s: V,
        ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>, tolerance: V) 
        -> Result<ExtendedMatrix<T, V>, String>
    {
        let c_matrix_multiplier_mem = young_modulus * thickness / (V::from(1f32) - poisson_ratio.my_powi(2));
        let mut c_matrix_mem = ExtendedMatrix::create(
            T::from(3u8), T::from(3u8), 
            vec![
                V::from(1f32), poisson_ratio, V::from(0f32),
                poisson_ratio, V::from(1f32), V::from(0f32),
                V::from(0f32), V::from(0f32), (V::from(1f32) - poisson_ratio) / V::from(2f32),
            ], 
            tolerance)?;
        c_matrix_mem.multiply_by_number(c_matrix_multiplier_mem);

        let rhs_matrix_mem = QuadFullPlateAuxFunctions::strain_displacement_matrix_mem(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, 
            ref_rotation_matrix, tolerance)?;
        let mut lhs_matrix_mem = rhs_matrix_mem.clone();
        lhs_matrix_mem.transpose();

        let mut matrix_mem = lhs_matrix_mem.multiply_by_matrix(&c_matrix_mem)?
            .multiply_by_matrix(&rhs_matrix_mem)?;
        matrix_mem.multiply_by_number(QuadFullPlateAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)? * alpha);
        

        let c_matrix_multiplier_bend = young_modulus * thickness.my_powi(3) / 
            (V::from(12f32) * (V::from(1f32) - poisson_ratio.my_powi(2)));
        let mut c_matrix_bend = ExtendedMatrix::create(
            T::from(3u8), T::from(3u8), 
            vec![
                V::from(1f32), poisson_ratio, V::from(0f32),
                poisson_ratio, V::from(1f32), V::from(0f32),
                V::from(0f32), V::from(0f32), (V::from(1f32) - poisson_ratio) / V::from(2f32),
            ], 
            tolerance)?;
        c_matrix_bend.multiply_by_number(c_matrix_multiplier_bend);

        let rhs_matrix_bend = QuadFullPlateAuxFunctions::strain_displacement_matrix_plate_bending(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, 
            ref_rotation_matrix, tolerance)?;
        let mut lhs_matrix_bend = rhs_matrix_bend.clone();
        lhs_matrix_bend.transpose();

        let mut matrix_bend = lhs_matrix_bend.multiply_by_matrix(&c_matrix_bend)?
            .multiply_by_matrix(&rhs_matrix_bend)?;
        matrix_bend.multiply_by_number(QuadFullPlateAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)? * alpha);

        let c_matrix_multiplier_shear = young_modulus * thickness * shear_factor / 
            (V::from(2f32) * (V::from(1f32) + poisson_ratio));
        let mut c_matrix_shear = ExtendedMatrix::create(
            T::from(2u8), T::from(2u8), 
            vec![
                V::from(1f32), V::from(0f32),
                V::from(0f32), V::from(1f32),
            ], 
            tolerance)?;
        c_matrix_shear.multiply_by_number(c_matrix_multiplier_shear * alpha);

        let rhs_matrix_shear = QuadFullPlateAuxFunctions::strain_displacement_matrix_plate_shear(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, 
            ref_rotation_matrix, tolerance)?;
        let mut lhs_matrix_shear = rhs_matrix_shear.clone();
        lhs_matrix_shear.transpose();

        let mut matrix_shear = lhs_matrix_shear.multiply_by_matrix(&c_matrix_shear)?
            .multiply_by_matrix(&rhs_matrix_shear)?;
        matrix_shear.multiply_by_number(QuadFullPlateAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, tolerance)? * alpha);

        let matrix = (matrix_mem.add_matrix(&matrix_bend)
            .map_err(|e| format!("QuadFullPlateAuxFunctions: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?).add_matrix(&matrix_shear)
            .map_err(|e| format!("QuadFullPlateAuxFunctions: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        Ok(matrix)
    }


    pub fn compose_node_dof_parameters<'a>(node_number: T)
        -> Result<Vec<DOFParameterData<T>>, &'a str>
    {
        let mut node_dof_parameters = Vec::new();
        for dof in 0..PLATE_NODE_DOF
        {
            let dof_parameter =
                GlobalDOFParameter::iterator().nth(dof)
                    .ok_or("QuadFullPlateAuxFunctions: Could not find dof parameter!")?;
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