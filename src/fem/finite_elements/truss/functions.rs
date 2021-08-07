use std::ops::{Add, Sub, Div, Rem, SubAssign, Mul, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::collections::HashMap;

use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::extract_element_value;

use crate::my_float::MyFloatTrait;

use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData, GlobalDOFParameter};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::truss::consts::{TRUSS2N2IP_NODES_NUMBER, TRUSS_NODE_DOF};
use crate::fem::finite_elements::functions::compare_with_tolerance;



pub struct TrussAuxFunctions<T, V>(T, V);


impl<T, V> TrussAuxFunctions<T, V>
    where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
             From<u8> + 'static,
          V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + From<f32> +
             Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
             MyFloatTrait + PartialOrd + 'static,
{
    fn length(node_1_number: T, node_2_number: T, nodes: &HashMap<T, FENode<V>>) -> V
    {
        let node_1 = nodes.get(&node_1_number).unwrap();
        let node_2 = nodes.get(&node_2_number).unwrap();

        ((node_1.x() - node_2.x()).my_powi(2) + (node_1.y() - node_2.y()).my_powi(2) +
        (node_1.z() - node_2.z()).my_powi(2)).my_sqrt()
    }


    pub fn nodes_number() -> T
    {
        let mut n = T::from(0u8);
        (0..TRUSS2N2IP_NODES_NUMBER).for_each(|_| n += T::from(1u8));
        n
    }


    pub fn node_dof() -> T
    {
        let mut n = T::from(0u8);
        (0..TRUSS_NODE_DOF).for_each(|_| n += T::from(1u8));
        n
    }


    pub fn rotation_matrix(node_1_number: T, node_2_number: T, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> ExtendedMatrix<T, V>
    {
        let node_1 = nodes.get(&node_1_number).unwrap();
        let node_2 = nodes.get(&node_2_number).unwrap();

        let x = node_2.x() - node_1.x();
        let y = node_2.y() - node_1.y();
        let z = node_2.z() - node_1.z();

        let length = TrussAuxFunctions::<T, V>::length(node_1_number, node_2_number, nodes);

        let (u, v, w) = (length, V::from(0f32), V::from(0f32));
        let alpha = ((x * u + y * v + z * w) / (length * length)).my_acos();
        let (rotation_axis_coord_x, mut rotation_axis_coord_y,
            mut rotation_axis_coord_z) = (V::from(0f32), V::from(0f32), V::from(0f32));
        if x != V::from(0f32) && y == V::from(0f32) && z == V::from(0f32)
        {
            rotation_axis_coord_z = x;
        }
        else
        {
            rotation_axis_coord_y = z * length;
            rotation_axis_coord_z = y * V::from(-1f32) * length;
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

        ExtendedMatrix::create(
            TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
            TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
            vec![
                [q_11, q_12, q_13], [V::from(0f32); TRUSS_NODE_DOF],
                [q_21, q_22, q_23], [V::from(0f32); TRUSS_NODE_DOF],
                [q_31, q_32, q_33], [V::from(0f32); TRUSS_NODE_DOF],
                [V::from(0f32); TRUSS_NODE_DOF], [q_11, q_12, q_13],
                [V::from(0f32); TRUSS_NODE_DOF], [q_21, q_22, q_23],
                [V::from(0f32); TRUSS_NODE_DOF], [q_31, q_32, q_33],
            ].concat(),
            tolerance,
        )
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


    fn dx_dr(x_1: V, x_2: V, r: V) -> V
    {
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.5f32), V::from(0f32),
            0) -
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.5f32), r, 1) +
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.5f32), V::from(0f32),
            0) +
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.5f32), r, 1)
    }


    fn jacobian(node_1_number: T, node_2_number: T, r: V, nodes: &HashMap<T, FENode<V>>) -> V
    {
        let length = TrussAuxFunctions::length(node_1_number, node_2_number, nodes);
        let x_1 = V::from(-1f32) * length / V::from(2f32);
        let x_2 = length / V::from(2f32);
        TrussAuxFunctions::<T, V>::dx_dr(x_1, x_2, r)
    }


    fn inverse_jacobian(node_1_number: T, node_2_number: T, r: V, nodes: &HashMap<T, FENode<V>>) -> V
    {
        V::from(1f32) / TrussAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    }


    fn determinant_of_jacobian(node_1_number: T, node_2_number: T, r: V,
        nodes: &HashMap<T, FENode<V>>) -> V
    {
        TrussAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    }


    fn dh1_dr(r: V) -> V
    {
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) -
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), r, 1)
    }


    fn dh2_dr(r: V) -> V
    {
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) +
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), r, 1)
    }


    pub fn strain_displacement_matrix(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> ExtendedMatrix<T, V>
    {
        let elements = vec![TrussAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32),
            V::from(0f32), TrussAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32), V::from(0f32)];
        let mut matrix = ExtendedMatrix::create(T::from(1u8),
            TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
            elements, tolerance);
        let inverse_jacobian = TrussAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        matrix.multiply_by_number(inverse_jacobian);
        matrix
    }


    pub fn area(area_1: V, area_2: Option<V>, r: V) -> V
    {
        if let Some(area_2) = area_2
        {
            (area_2 - area_1) / V::from(2f32) * r + area_1 -
                (area_2 - area_1) / V::from(2f32) * V::from(-1f32)
        }
        else
        {
            area_1
        }
    }


    pub fn local_stiffness_matrix(node_1_number: T, node_2_number: T, young_modulus: V, area_1: V,
        area_2: Option<V>, alpha: V, r: V, local_stiffness_matrix: &ExtendedMatrix<T, V>,
        tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let current_area = TrussAuxFunctions::<T, V>::area(area_1, area_2, r);

        let mut lhs_matrix = TrussAuxFunctions::strain_displacement_matrix(
            node_1_number, node_2_number, r, tolerance, nodes);

        lhs_matrix.transpose();

        lhs_matrix.multiply_by_number(young_modulus * current_area);

        let rhs_matrix = TrussAuxFunctions::strain_displacement_matrix(
                node_1_number, node_2_number, r, tolerance, nodes);

        return match lhs_matrix.multiply_by_matrix(&rhs_matrix)
        {
            Ok(mut matrix) =>
                {
                    matrix.multiply_by_number(TrussAuxFunctions::determinant_of_jacobian(
                        node_1_number, node_2_number, r, nodes) * alpha);

                    match local_stiffness_matrix.add_matrix(&matrix)
                    {
                        Ok(matrix) => Ok(matrix),
                        Err(e) =>
                            Err(format!("Truss2n2ip: Local stiffness matrix cannot be \
                                calculated! Reason: {}", e)),
                    }
                },
            Err(e) => Err(format!("Truss2n2ip: Local stiffness matrix cannot be \
                                calculated! Reason: {}", e)),
        }
    }


    pub fn compose_node_dof_parameters<'a>(node_number: T)
        -> Result<Vec<DOFParameterData<T>>, &'a str>
    {
        let mut node_dof_parameters = Vec::new();
        for dof in 0..TRUSS_NODE_DOF
        {
            let dof_parameter =
                GlobalDOFParameter::iterator().nth(dof)
                    .ok_or("Truss2n2ip: Could not find dof parameter!")?;
            let dof_parameter = DOFParameterData { node_number,
                dof_parameter: *dof_parameter
            };
            node_dof_parameters.push(dof_parameter);
        }
        Ok(node_dof_parameters)
    }


    pub fn extract_column_matrix_values(column_matrix: &ExtendedMatrix<T, V>) -> Vec<V>
    {
        let mut values = Vec::new();
        let shape = column_matrix.get_shape();
        let all_values = column_matrix.extract_all_elements_values();

        let mut row = T::from(0u8);
        while row < shape.0
        {
            let mut column = T::from(0u8);
            while column < shape.1
            {
                let value = extract_element_value(row, column, &all_values);
                values.push(value);
                column += T::from(1u8);
            }
            row += T::from(1u8);
        }
        values
    }
}
