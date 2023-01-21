use std::collections::HashMap;

use extended_matrix::
{
    Matrix, FloatTrait, Vector3, VectorTrait, Position, BasicOperationsTrait, SquareMatrixTrait, SquareMatrix,
};

use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData, GlobalDOFParameter};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::membrane::consts::{MEM4N4IP_NODES_NUMBER, MEMBRANE_NODE_DOF};
use crate::fem::finite_elements::functions::extract_unique_elements_of_rotation_matrix;


pub struct QuadFullMemAuxFunctions<V>(V);


impl<V> QuadFullMemAuxFunctions<V>
    where V: FloatTrait<Output = V>
{
    pub fn nodes_number() -> usize
    {
        MEM4N4IP_NODES_NUMBER
    }


    pub fn node_dof() -> usize
    {
        MEMBRANE_NODE_DOF
    }


    pub fn rotation_matrix(
        node_2_number: u32, 
        node_3_number: u32, 
        node_4_number: u32, 
        ref_nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    ) 
        -> Result<Matrix<V>, String>
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

        let edge_3_4_vector = Vector3::create(&[edge_3_4_x, edge_3_4_y, edge_3_4_z]);
        let edge_3_2_vector = Vector3::create(&[edge_3_2_x, edge_3_2_y, edge_3_2_z]);

        let normal_through_node_3 = edge_3_4_vector.cross_product(&edge_3_2_vector);
        let normal_through_node_3_length = normal_through_node_3.norm()?;

        let direction_vector = Vector3::create(&[0.0, 0.0, normal_through_node_3_length]);

        let shrinked_rotation_matrix = normal_through_node_3
            .rotation_matrix_to_align_with_vector(&direction_vector, rel_tol, abs_tol)?;

        let q_11 = shrinked_rotation_matrix.get_element_value(&Position(0, 0))?;
        let q_12 = shrinked_rotation_matrix.get_element_value(&Position(0, 1))?;
        let q_13 = shrinked_rotation_matrix.get_element_value(&Position(0, 2))?;
        let q_21 = shrinked_rotation_matrix.get_element_value(&Position(1, 0))?;
        let q_22 = shrinked_rotation_matrix.get_element_value(&Position(1, 1))?;
        let q_23 = shrinked_rotation_matrix.get_element_value(&Position(1, 2))?;
        let q_31 = shrinked_rotation_matrix.get_element_value(&Position(2, 0))?;
        let q_32 = shrinked_rotation_matrix.get_element_value(&Position(2, 1))?;
        let q_33 = shrinked_rotation_matrix.get_element_value(&Position(2, 2))?;

        // let normal_through_node_3_x = edge_3_4_y * edge_3_2_z - edge_3_4_z * edge_3_2_y;
        // let normal_through_node_3_y = edge_3_4_z * edge_3_2_x - edge_3_4_x * edge_3_2_z;
        // let normal_through_node_3_z = edge_3_4_x * edge_3_2_y - edge_3_4_y * edge_3_2_x;

        // let normal_through_node_3_length = (
        //         (normal_through_node_3_x).my_powi(2) + 
        //         (normal_through_node_3_y).my_powi(2) + 
        //         (normal_through_node_3_z).my_powi(2)
        //     )
        //     .my_sqrt();

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
        
        let rotation_matrix = Matrix::create(
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
            &[
                [q_11, q_12, q_13], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [q_21, q_22, q_23], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [q_31, q_32, q_33], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],

                [&V::from(0f32); MEMBRANE_NODE_DOF], [q_11, q_12, q_13],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [q_21, q_22, q_23],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [q_31, q_32, q_33],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],

                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [q_11, q_12, q_13], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [q_21, q_22, q_23], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [q_31, q_32, q_33], [&V::from(0f32); MEMBRANE_NODE_DOF],

                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [q_11, q_12, q_13],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [q_21, q_22, q_23],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [&V::from(0f32); MEMBRANE_NODE_DOF],
                [&V::from(0f32); MEMBRANE_NODE_DOF], [q_31, q_32, q_33],
            ]
            .concat(),
        );
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
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_1 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_1 * V::from(0.25f32), r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_1 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_2 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_2 * V::from(0.25f32), r, 1) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_2 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_3 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_3 * V::from(0.25f32), r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_3 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_4 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_4 * V::from(0.25f32), r, 1) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_4 * V::from(0.25f32) * s, r, 1)
    }


    fn dx_ds(x_1: V, x_2: V, x_3: V, x_4: V, r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_1 * V::from(0.25f32), s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_1 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_1 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_2 * V::from(0.25f32), s, 1) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_2 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_2 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_3 * V::from(0.25f32), s, 1) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_3 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_3 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_4 * V::from(0.25f32), s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_4 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, x_4 * V::from(0.25f32) * r, s, 1)
    }


    fn dy_dr(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_1 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_1 * V::from(0.25f32), r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_1 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_2 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_2 * V::from(0.25f32), r, 1) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_2 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_3 * V::from(0.25f32) * s, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_3 * V::from(0.25f32), r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_3 * V::from(0.25f32) * s, r, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_4 * V::from(0.25f32) * s, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_4 * V::from(0.25f32), r, 1) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_4 * V::from(0.25f32) * s, r, 1)
    }


    fn dy_ds(y_1: V, y_2: V, y_3: V, y_4: V, r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_1 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_1 * V::from(0.25f32), s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_1 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_1 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_2 * V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_2 * V::from(0.25f32), s, 1) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_2 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_2 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_3 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_3 * V::from(0.25f32), s, 1) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_3 * V::from(0.25f32) * r, V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_3 * V::from(0.25f32) * r, s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_4 * V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_4 * V::from(0.25f32), s, 1) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_4 * V::from(0.25f32) * r, V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, y_4 * V::from(0.25f32) * r, s, 1)
    }


    // fn jacobian(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T,
    //     r: V, s: V, ref_nodes: &HashMap<T, FENode<V>>, ref_rotation_matrix: &ExtendedMatrix<T, V>,
    //     tolerance: V) -> Result<ExtendedMatrix<T, V>, String>
    pub fn jacobian(
        node_1_number: u32, 
        node_2_number: u32, 
        node_3_number: u32, 
        node_4_number: u32,
        r: V, 
        s: V, 
        ref_nodes: &HashMap<u32, FENode<V>>, 
        ref_rotation_matrix: &Matrix<V>,
    ) 
        -> Result<SquareMatrix<V>, String>
    {
        let unique_elements_of_rotation_matrix = extract_unique_elements_of_rotation_matrix::<V>(
            ref_rotation_matrix,
        )?;
        let shrinked_rotation_matrix = Matrix::create(
            3,
            3, 
            &unique_elements_of_rotation_matrix,
        );

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
        
        let node_1_direction = Matrix::create(
            3,
            1,
            &node_1_direction_vector,
        );
        let node_2_direction = Matrix::create(
            3,
            1,
            &node_2_direction_vector,
        );
        let node_4_direction = Matrix::create(
            3,
            1,
            &node_4_direction_vector,
        );

        let transformed_node_1_direction = shrinked_rotation_matrix.multiply(&node_1_direction)?;
        let transformed_node_2_direction = shrinked_rotation_matrix.multiply(&node_2_direction)?;
        let transformed_node_4_direction = shrinked_rotation_matrix.multiply(&node_4_direction)?;

        let transformed_node_1_direction_x = transformed_node_1_direction.get_element_value(&Position(0, 0))?;
        let transformed_node_1_direction_y = transformed_node_1_direction.get_element_value(&Position(1, 0))?;
        let transformed_node_1_direction_z = transformed_node_1_direction.get_element_value(&Position(2, 0))?;
        let transformed_node_2_direction_x = transformed_node_2_direction.get_element_value(&Position(0, 0))?;
        let transformed_node_2_direction_y = transformed_node_2_direction.get_element_value(&Position(1, 0))?;
        let transformed_node_2_direction_z = transformed_node_2_direction.get_element_value(&Position(2, 0))?;
        let transformed_node_4_direction_x = transformed_node_4_direction.get_element_value(&Position(0, 0))?;
        let transformed_node_4_direction_y = transformed_node_4_direction.get_element_value(&Position(1, 0))?;
        let transformed_node_4_direction_z = transformed_node_4_direction.get_element_value(&Position(2, 0))?;

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
            QuadFullMemAuxFunctions::<V>::dx_dr(
                transformed_node_1_direction_x, 
                transformed_node_2_direction_x, 
                V::from(0f32), 
                transformed_node_4_direction_x, 
                r,
                s,
            ),
            QuadFullMemAuxFunctions::<V>::dy_dr(
                transformed_node_1_direction_y,
                transformed_node_2_direction_y, 
                V::from(0f32),
                transformed_node_4_direction_y,
                r,
                s,
            ),
            QuadFullMemAuxFunctions::<V>::dx_ds(
                transformed_node_1_direction_x,
                transformed_node_2_direction_x, 
                V::from(0f32),
                transformed_node_4_direction_x,
                r,
                s,
            ),
            QuadFullMemAuxFunctions::<V>::dy_ds(
                transformed_node_1_direction_y,
                transformed_node_2_direction_y, 
                V::from(0f32),
                transformed_node_4_direction_y,
                r,
                s,
            ),
        ];

        let jacobian = SquareMatrix::create(2, &jacobian_elements);

        Ok(jacobian)
    }


    fn inverse_jacobian(
        node_1_number: u32, 
        node_2_number: u32, 
        node_3_number: u32,
        node_4_number: u32,
        r: V, 
        s: V, 
        ref_nodes: &HashMap<u32, FENode<V>>, 
        ref_rotation_matrix: &Matrix<V>,
        rel_tol: V,
    ) 
        -> Result<Matrix<V>, String>
    {
        let jacobian = QuadFullMemAuxFunctions::<V>::jacobian(
            node_1_number, 
            node_2_number, 
            node_3_number, 
            node_4_number, 
            r, 
            s, 
            ref_nodes, 
            ref_rotation_matrix,
        )?;
        let mut x = vec![V::from(0f32); jacobian.get_shape().0];
        let inverse_jacobian = jacobian.inverse(&mut x, rel_tol)?;
        Ok(inverse_jacobian)
    }


    fn determinant_of_jacobian(
        node_1_number: u32,
        node_2_number: u32,
        node_3_number: u32,
        node_4_number: u32,
        r: V, 
        s: V,
        ref_nodes: &HashMap<u32, FENode<V>>,
        ref_rotation_matrix: &Matrix<V>,
        rel_tol: V,
    ) 
        -> Result<V, String>
    {
        let jacobian = QuadFullMemAuxFunctions::<V>::jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix
        )?;
        let determinant_of_jacobian = jacobian.determinant(rel_tol);
        Ok(determinant_of_jacobian)
    }


    fn dh1_dr(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) + 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), r, 1) + 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh2_dr(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) - 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), r, 1) - 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh3_dr(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) - 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), r, 1) + 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh4_dr(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * s, V::from(0f32), 0) + 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), r, 1) - 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * s, r, 1)
    }


    fn dh1_ds(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), s, 1) + 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) + 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh2_ds(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) +
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), s, 1) - 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) - 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh3_ds(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), s, 1) - 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) + 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh4_ds(r: V, s: V) -> V
    {
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), V::from(0f32), 0) -
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32), s, 1) + 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * r, V::from(0f32), 0) - 
        QuadFullMemAuxFunctions::<V>::derivative_x(
            QuadFullMemAuxFunctions::<V>::power_func_x, V::from(0.25f32) * r, s, 1)
    }


    fn dh_dx_dh_dy(
        node_1_number: u32, 
        node_2_number: u32, 
        node_3_number: u32, 
        node_4_number: u32,
        r: V, 
        s: V, 
        ref_nodes: &HashMap<u32, FENode<V>>, 
        ref_rotation_matrix: &Matrix<V>,
        rel_tol: V,
    ) 
        -> Result<Matrix<V>, String>
    {
        let inverse_jacobian = QuadFullMemAuxFunctions::inverse_jacobian(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, rel_tol,
        )?;
        let dh_dr_dh_ds = Matrix::create(
            2,
            4, 
            &[
                QuadFullMemAuxFunctions::<V>::dh1_dr(r, s), QuadFullMemAuxFunctions::<V>::dh2_dr(r, s), 
                QuadFullMemAuxFunctions::<V>::dh3_dr(r, s), QuadFullMemAuxFunctions::<V>::dh4_dr(r, s),
                QuadFullMemAuxFunctions::<V>::dh1_ds(r, s), QuadFullMemAuxFunctions::<V>::dh2_ds(r, s), 
                QuadFullMemAuxFunctions::<V>::dh3_ds(r, s), QuadFullMemAuxFunctions::<V>::dh4_ds(r, s),
            ], 
        );
        Ok(inverse_jacobian.multiply_by_matrix(&dh_dr_dh_ds)?)

    }


    pub fn strain_displacement_matrix(
        node_1_number: u32, 
        node_2_number: u32, 
        node_3_number: u32, 
        node_4_number: u32,
        r: V, 
        s: V, 
        ref_nodes: &HashMap<u32, FENode<V>>, 
        ref_rotation_matrix: &Matrix<V>,
        rel_tol: V,
    ) 
        -> Result<Matrix<V>, String>
    {
        let dh_dx_dh_dy_matrix = QuadFullMemAuxFunctions::<V>::dh_dx_dh_dy(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, rel_tol,
        )?;

        let dh1_dx = dh_dx_dh_dy_matrix.get_element_value(&Position(0, 0))?;
        let dh2_dx = dh_dx_dh_dy_matrix.get_element_value(&Position(0, 1))?;
        let dh3_dx = dh_dx_dh_dy_matrix.get_element_value(&Position(0, 2))?;
        let dh4_dx = dh_dx_dh_dy_matrix.get_element_value(&Position(0, 3))?;

        let dh1_dy = dh_dx_dh_dy_matrix.get_element_value(&Position(1, 0))?;
        let dh2_dy = dh_dx_dh_dy_matrix.get_element_value(&Position(1, 1))?;
        let dh3_dy = dh_dx_dh_dy_matrix.get_element_value(&Position(1, 2))?;
        let dh4_dy = dh_dx_dh_dy_matrix.get_element_value(&Position(1, 3))?;

        let elements = [
            *dh1_dx, V::from(0f32), V::from(0f32), *dh2_dx, V::from(0f32), V::from(0f32), 
            *dh3_dx, V::from(0f32), V::from(0f32), *dh4_dx, V::from(0f32), V::from(0f32), 

            V::from(0f32), *dh1_dy, V::from(0f32), V::from(0f32), *dh2_dy, V::from(0f32),
            V::from(0f32), *dh3_dy, V::from(0f32), V::from(0f32), *dh4_dy, V::from(0f32),

            *dh1_dy, *dh1_dx, V::from(0f32), *dh2_dy, *dh2_dx, V::from(0f32), 
            *dh3_dy, *dh3_dx, V::from(0f32), *dh4_dy, *dh4_dx, V::from(0f32),
        ];

        let matrix = Matrix::create(
            3, 
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(), 
            &elements,
        );

        Ok(matrix)
    }


    pub fn local_stiffness_matrix(
        node_1_number: u32, 
        node_2_number: u32, 
        node_3_number: u32, 
        node_4_number: u32,
        young_modulus: V, 
        poisson_ratio: V, 
        thickness: V, 
        alpha: V, 
        r: V, 
        s: V, 
        ref_local_stiffness_matrix: &Matrix<V>, 
        ref_nodes: &HashMap<u32, FENode<V>>, 
        ref_rotation_matrix: &Matrix<V>, 
        rel_tol: V,
    )
        -> Result<Matrix<V>, String>
    {
        let c_matrix_multiplier = young_modulus * thickness / (V::from(1f32) - poisson_ratio.my_powi(2));
        let mut c_matrix = Matrix::create(
                3,
                3, 
                &[
                    V::from(1f32), poisson_ratio, V::from(0f32),
                    poisson_ratio, V::from(1f32), V::from(0f32),
                    V::from(0f32), V::from(0f32), (V::from(1f32) - poisson_ratio) / V::from(2f32),
                ],     
            )
            .multiply_by_scalar(c_matrix_multiplier);

        let mut lhs_matrix = QuadFullMemAuxFunctions::strain_displacement_matrix(
                node_1_number, 
                node_2_number, 
                node_3_number, 
                node_4_number, 
                r, 
                s, 
                ref_nodes, 
                ref_rotation_matrix, 
                rel_tol,
            )?
            .transpose();

        let rhs_matrix = QuadFullMemAuxFunctions::strain_displacement_matrix(
            node_1_number, node_2_number, node_3_number, node_4_number, r, s, ref_nodes, ref_rotation_matrix, rel_tol,
        )?;

        return match (lhs_matrix.multiply(&c_matrix)?).multiply(&rhs_matrix)
        {
            Ok(mut matrix) =>
                {
                    matrix.multiply_by_number(QuadFullMemAuxFunctions::determinant_of_jacobian(
                            node_1_number, 
                            node_2_number, 
                            node_3_number, 
                            node_4_number, 
                            r, 
                            s, 
                            ref_nodes, 
                            ref_rotation_matrix, 
                            rel_tol,
                        )? * 
                        alpha,
                    );

                    match ref_local_stiffness_matrix.add(&matrix)
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


    pub fn compose_node_dof_parameters(node_number: u32) -> Result<Vec<DOFParameterData>, String>
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


    pub fn extract_column_matrix_values(column_matrix: &Matrix<V>) -> Result<Vec<V>, String>
    {
        let mut values = Vec::new();
        let shape = column_matrix.copy_shape();

        let mut row = 0;
        while row < shape.0
        {
            let mut column = 0;
            while column < shape.1
            {
                let value = column_matrix.get_element_value(&Position(row, column))?;
                values.push(value);
                column += 1;
            }
            row += 1;
        }
        Ok(values)
    }
}
