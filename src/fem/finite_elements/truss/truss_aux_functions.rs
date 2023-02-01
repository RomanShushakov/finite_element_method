use std::collections::HashMap;

use extended_matrix::{Matrix, FloatTrait, Vector3, VectorTrait, BasicOperationsTrait, Position};

use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData, GlobalDOFParameter};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::truss::consts::{TRUSS2N2IP_NODES_NUMBER, TRUSS_NODE_DOF};
// use crate::fem::finite_elements::functions::compare_with_tolerance;



pub struct TrussAuxFunctions<V>(V);


impl<V> TrussAuxFunctions<V>
    where V: FloatTrait<Output = V>
{
    fn length(
        node_1_number: u32,
        node_2_number: u32,
        nodes: &HashMap<u32, FENode<V>>
    ) 
        -> V
    {
        let node_1 = nodes.get(&node_1_number).unwrap();
        let node_2 = nodes.get(&node_2_number).unwrap();

        (
            (node_1.copy_x() - node_2.copy_x()).my_powi(2) + 
            (node_1.copy_y() - node_2.copy_y()).my_powi(2) +
            (node_1.copy_z() - node_2.copy_z()).my_powi(2)
        ).my_sqrt()
    }


    pub fn nodes_number() -> usize
    {
        TRUSS2N2IP_NODES_NUMBER
    }


    pub fn node_dof() -> usize
    {
        TRUSS_NODE_DOF
    }


    pub fn rotation_matrix(
        node_1_number: u32,
        node_2_number: u32,
        nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    ) 
        -> Result<Matrix<V>, String>
    {
        let node_1 = nodes.get(&node_1_number).unwrap();
        let node_2 = nodes.get(&node_2_number).unwrap();

        let x = node_2.copy_x() - node_1.copy_x();
        let y = node_2.copy_y() - node_1.copy_y();
        let z = node_2.copy_z() - node_1.copy_z();

        let truss_element_vector = Vector3::create(&[x, y, z]);
        let truss_element_length = truss_element_vector.norm()?;
        let direction_vector = Vector3::create(
            &[truss_element_length, V::from(0f32), V::from(0f32)],
        );

        let shrinked_rotation_matrix = truss_element_vector
            .rotation_matrix_to_align_with_vector(&direction_vector, rel_tol, abs_tol)?;

        let q_11 = *shrinked_rotation_matrix.get_element_value(&Position(0, 0))?;
        let q_12 = *shrinked_rotation_matrix.get_element_value(&Position(0, 1))?;
        let q_13 = *shrinked_rotation_matrix.get_element_value(&Position(0, 2))?;
        let q_21 = *shrinked_rotation_matrix.get_element_value(&Position(1, 0))?;
        let q_22 = *shrinked_rotation_matrix.get_element_value(&Position(1, 1))?;
        let q_23 = *shrinked_rotation_matrix.get_element_value(&Position(1, 2))?;
        let q_31 = *shrinked_rotation_matrix.get_element_value(&Position(2, 0))?;
        let q_32 = *shrinked_rotation_matrix.get_element_value(&Position(2, 1))?;
        let q_33 = *shrinked_rotation_matrix.get_element_value(&Position(2, 2))?;

        let rotation_matrix = Matrix::create(
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            &vec![
                [q_11, q_12, q_13], [V::from(0f32); TRUSS_NODE_DOF],
                [q_21, q_22, q_23], [V::from(0f32); TRUSS_NODE_DOF],
                [q_31, q_32, q_33], [V::from(0f32); TRUSS_NODE_DOF],
                [V::from(0f32); TRUSS_NODE_DOF], [q_11, q_12, q_13],
                [V::from(0f32); TRUSS_NODE_DOF], [q_21, q_22, q_23],
                [V::from(0f32); TRUSS_NODE_DOF], [q_31, q_32, q_33],
            ].concat(),
        );

        // let length = TrussAuxFunctions::<V>::length(node_1_number, node_2_number, nodes);

        // let (u, v, w) = (length, V::from(0f32), V::from(0f32));
        // let alpha = ((x * u + y * v + z * w) / (length * length)).my_acos();
        // let (rotation_axis_coord_x, mut rotation_axis_coord_y,
        //     mut rotation_axis_coord_z) = (V::from(0f32), V::from(0f32), V::from(0f32));
        // if x != V::from(0f32) && y == V::from(0f32) && z == V::from(0f32)
        // {
        //     rotation_axis_coord_z = x;
        // }
        // else
        // {
        //     rotation_axis_coord_y = z * length;
        //     rotation_axis_coord_z = y * V::from(-1f32) * length;
        // }
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

        // let rotation_matrix = ExtendedMatrix::create(
        //     TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
        //     TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
        //     vec![
        //         [q_11, q_12, q_13], [V::from(0f32); TRUSS_NODE_DOF],
        //         [q_21, q_22, q_23], [V::from(0f32); TRUSS_NODE_DOF],
        //         [q_31, q_32, q_33], [V::from(0f32); TRUSS_NODE_DOF],
        //         [V::from(0f32); TRUSS_NODE_DOF], [q_11, q_12, q_13],
        //         [V::from(0f32); TRUSS_NODE_DOF], [q_21, q_22, q_23],
        //         [V::from(0f32); TRUSS_NODE_DOF], [q_31, q_32, q_33],
        //     ].concat(),
        //     tolerance,
        // )?;

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


    fn dx_dr(x_1: V, x_2: V, r: V) -> V
    {
        TrussAuxFunctions::<V>::derivative_x(
            TrussAuxFunctions::<V>::power_func_x, x_1 * V::from(0.5f32), V::from(0f32),
            0) -
        TrussAuxFunctions::<V>::derivative_x(
            TrussAuxFunctions::<V>::power_func_x, x_1 * V::from(0.5f32), r, 1) +
        TrussAuxFunctions::<V>::derivative_x(
            TrussAuxFunctions::<V>::power_func_x, x_2 * V::from(0.5f32), V::from(0f32),
            0) +
        TrussAuxFunctions::<V>::derivative_x(
            TrussAuxFunctions::<V>::power_func_x, x_2 * V::from(0.5f32), r, 1)
    }


    fn jacobian(node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, FENode<V>>) -> V
    {
        let length = TrussAuxFunctions::length(node_1_number, node_2_number, nodes);
        let x_1 = V::from(-1f32) * length / V::from(2f32);
        let x_2 = length / V::from(2f32);
        TrussAuxFunctions::<V>::dx_dr(x_1, x_2, r)
    }


    fn inverse_jacobian(node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, FENode<V>>)
        -> V
    {
        V::from(1f32) / TrussAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    }


    fn determinant_of_jacobian(node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, FENode<V>>) -> V
    {
        TrussAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    }


    fn dh1_dr(r: V) -> V
    {
        TrussAuxFunctions::<V>::derivative_x(
            TrussAuxFunctions::<V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) -
        TrussAuxFunctions::<V>::derivative_x(
            TrussAuxFunctions::<V>::power_func_x, V::from(0.5f32), r, 1)
    }


    fn dh2_dr(r: V) -> V
    {
        TrussAuxFunctions::<V>::derivative_x(
            TrussAuxFunctions::<V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) +
        TrussAuxFunctions::<V>::derivative_x(
            TrussAuxFunctions::<V>::power_func_x, V::from(0.5f32), r, 1)
    }


    pub fn strain_displacement_matrix(
        node_1_number: u32, 
        node_2_number: u32, 
        r: V, 
        nodes: &HashMap<u32, FENode<V>>
    ) 
        -> Result<Matrix<V>, String>
    {
        let elements = vec![
            TrussAuxFunctions::<V>::dh1_dr(r), V::from(0f32), V::from(0f32), 
            TrussAuxFunctions::<V>::dh2_dr(r), V::from(0f32), V::from(0f32),
        ];
        let mut matrix = Matrix::create(
            1,
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            &elements,
        );
        let inverse_jacobian = TrussAuxFunctions::inverse_jacobian(node_1_number, node_2_number, r, nodes);
        matrix.multiply_by_scalar(inverse_jacobian);
        Ok(matrix)
    }


    pub fn area(area_1: V, area_2: Option<V>, r: V) -> V
    {
        if let Some(area_2) = area_2
        {
            (area_2 - area_1) / V::from(2f32) * r + area_1 - (area_2 - area_1) / V::from(2f32) * V::from(-1f32)
        }
        else
        {
            area_1
        }
    }


    pub fn local_stiffness_matrix(
        node_1_number: u32, 
        node_2_number: u32, 
        young_modulus: V, 
        area_1: V,
        area_2: Option<V>, 
        alpha: V, 
        r: V, 
        local_stiffness_matrix: &Matrix<V>,
        nodes: &HashMap<u32, FENode<V>>,
    ) 
        -> Result<Matrix<V>, String>
    {
        let current_area = TrussAuxFunctions::<V>::area(area_1, area_2, r);

        let mut lhs_matrix = TrussAuxFunctions::strain_displacement_matrix(
            node_1_number, node_2_number, r, nodes,
        )?;

        lhs_matrix.transpose();

        lhs_matrix.multiply_by_scalar(young_modulus * current_area);

        let rhs_matrix = TrussAuxFunctions::strain_displacement_matrix(
            node_1_number, node_2_number, r, nodes,
        )?;

        return match lhs_matrix.multiply(&rhs_matrix)
        {
            Ok(mut matrix) =>
                {
                    matrix.multiply_by_scalar(
                        TrussAuxFunctions::determinant_of_jacobian(node_1_number, node_2_number, r, nodes) * alpha,
                    );

                    match local_stiffness_matrix.add(&matrix)
                    {
                        Ok(matrix) => Ok(matrix),
                        Err(e) => Err(format!("Truss element: Local stiffness matrix cannot be calculated! Reason: {e}")),
                    }
                },
            Err(e) => Err(format!("Truss element: Local stiffness matrix cannot be calculated! Reason: {e}")),
        }
    }


    pub fn compose_node_dof_parameters(node_number: u32) -> Result<Vec<DOFParameterData>, String>
    {
        let mut node_dof_parameters = Vec::new();
        for dof in 0..TRUSS_NODE_DOF
        {
            let dof_parameter = GlobalDOFParameter::iterator()
                .nth(dof)
                .ok_or("Truss2n2ip: Could not find dof parameter!".to_string())?;
            let dof_parameter = DOFParameterData::create(node_number, *dof_parameter);
            node_dof_parameters.push(dof_parameter);
        }
        Ok(node_dof_parameters)
    }


    pub fn extract_column_matrix_values(column_matrix: &Matrix<V>) -> Result<Vec<V>, String>
    {
        let mut values = Vec::new();
        let shape = column_matrix.get_shape();

        let mut row = 0;
        while row < shape.0
        {
            let mut column = 0;
            while column < shape.1
            {
                let value =  column_matrix.get_element_value(&Position(row, column))?;
                values.push(*value);
                column += 1;
            }
            row += 1;
        }
        Ok(values)
    }
}
