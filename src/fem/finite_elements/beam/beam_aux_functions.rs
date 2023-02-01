use std::collections::HashMap;
use std::f32::consts::PI;

use extended_matrix::{Matrix, FloatTrait, Vector3, Position, BasicOperationsTrait, VectorTrait};

use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData, GlobalDOFParameter};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::beam::consts::{BEAM2N1IPT_NODES_NUMBER, BEAM_NODE_DOF};
use crate::fem::finite_elements::functions::compare_with_tolerance;


pub struct BeamAuxFunctions<V>(V);


impl<V> BeamAuxFunctions<V>
    where V: FloatTrait<Output = V>
{
    pub fn length(node_1_number: u32, node_2_number: u32, nodes: &HashMap<u32, FENode<V>>) -> V
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
        BEAM2N1IPT_NODES_NUMBER
    }


    pub fn node_dof() -> usize
    {
        BEAM_NODE_DOF
    }


    fn find_components_of_line_a_perpendicular_to_line_b(
        line_a: &[V; 3], 
        line_b: &[V; 3], 
    )
        -> Result<[V; 3], String>
    {
        let vector_a = Vector3::create(line_a);
        let vector_b = Vector3::create(line_b);
        let projection = vector_a.projection_perpendicular_to_vector(&vector_b);

        // let a_x = V::from(-1f32) * line_a[0];
        // let a_y = V::from(-1f32) * line_a[1];
        // let a_z = V::from(-1f32) * line_a[2];

        // let a = ExtendedMatrix::create(T::from(3u8),
        //     T::from(1u8), vec![a_x, a_y, a_z], tolerance)?;

        // let b_x = line_b[0];
        // let b_y = line_b[1];
        // let b_z = line_b[2];

        // let norm = V::from(1f32) / (b_x.my_powi(2) + b_y.my_powi(2) + b_z.my_powi(2));

        // let mut coeff_matrix = ExtendedMatrix::create(T::from(3u8),
        //     T::from(3u8), vec![
        //         V::from(-1f32) * b_z * b_z - b_y * b_y, b_x * b_y, b_x * b_z,
        //         b_y * b_x, V::from(-1f32) * b_x * b_x - b_z * b_z,	b_y * b_z,
        //         b_z * b_x,	b_z * b_y, V::from(-1f32) * b_y * b_y - b_x * b_x,
        //     ], tolerance)?;

        // coeff_matrix.multiply_by_number(norm);

        // let components_of_line_a_perpendicular_to_line_b_matrix = coeff_matrix
        //     .multiply_by_matrix(&a)?;

        // let a_perpendicular_to_b_x = matrix_element_value_extractor(
        //     T::from(0u8), T::from(0u8), &components_of_line_a_perpendicular_to_line_b_matrix)?;

        // let a_perpendicular_to_b_y = matrix_element_value_extractor(
        //     T::from(1u8), T::from(0u8), &components_of_line_a_perpendicular_to_line_b_matrix)?;

        // let a_perpendicular_to_b_z = matrix_element_value_extractor(
        //     T::from(2u8), T::from(0u8), &components_of_line_a_perpendicular_to_line_b_matrix)?;

        // Ok([a_perpendicular_to_b_x, a_perpendicular_to_b_y, a_perpendicular_to_b_z])

        Ok(projection.get_components())
    }


    pub fn rotation_matrix(
        node_1_number: u32, 
        node_2_number: u32, 
        local_axis_1_direction: &[V; 3],
        angle: V, 
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
        
        let beam_element_vector = Vector3::create(&[x, y, z]);
        let length = beam_element_vector.norm()?;
        let direction_vector = Vector3::create(&[length, V::from(0f32), V::from(0f32)]);

        let shrinked_rotation_matrix = beam_element_vector
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

        // let length = BeamAuxFunctions::<T, V>::length(node_1_number, node_2_number, nodes);

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
        //     rotation_axis_coord_z = V::from(-1f32) * y * length;
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

        let components_projection_of_beam_section_orientation_vector = 
            BeamAuxFunctions::<V>::find_components_of_line_a_perpendicular_to_line_b(
                local_axis_1_direction, &[x, y, z],
            )?;

        let interim_rotation_matrix = Matrix::create(
            3, 3, &[q_11, q_12, q_13, q_21, q_22, q_23, q_31, q_32, q_33],
        );

        let projection_of_beam_section_orientation = Matrix::create(
            3, 
            1,
            &components_projection_of_beam_section_orientation_vector,
        );

        let transformed_projection_of_beam_section_orientation = interim_rotation_matrix
            .multiply(&projection_of_beam_section_orientation)?;

        let transformed_projection_of_beam_section_orientation_x = 
            *transformed_projection_of_beam_section_orientation.get_element_value(&Position(0, 0))?;

        let transformed_projection_of_beam_section_orientation_y =
            *transformed_projection_of_beam_section_orientation.get_element_value(&Position(1, 0))?;

        let transformed_projection_of_beam_section_orientation_z =
            *transformed_projection_of_beam_section_orientation.get_element_value(&Position(2, 0))?;

        let angle_between_beam_section_local_axis_1_direction_and_axis_t =
            (transformed_projection_of_beam_section_orientation_z /
            (transformed_projection_of_beam_section_orientation_x.my_powi(2) +
            transformed_projection_of_beam_section_orientation_y.my_powi(2) +
            transformed_projection_of_beam_section_orientation_z.my_powi(2))
                .my_sqrt()
            ).my_acos();

        let total_angle = angle +
            angle_between_beam_section_local_axis_1_direction_and_axis_t;

        let c_x = compare_with_tolerance(x / length, abs_tol);
        let c_y = compare_with_tolerance(y / length, abs_tol);
        let c_z = compare_with_tolerance(z / length, abs_tol);
        let c_xz = compare_with_tolerance((c_x.my_powi(2) + c_z.my_powi(2)).my_sqrt(), abs_tol);

        let c = compare_with_tolerance(total_angle.my_cos(), abs_tol);
        let s = compare_with_tolerance(total_angle.my_sin(), abs_tol);

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

        let rotation_matrix = Matrix::create(
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &[
                [r_11, r_12, r_13], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [r_21, r_22, r_23], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [r_31, r_32, r_33], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],

                [V::from(0f32); BEAM_NODE_DOF / 2], [r_11, r_12, r_13],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [r_21, r_22, r_23],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [r_31, r_32, r_33],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],

                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [r_11, r_12, r_13], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [r_21, r_22, r_23], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [r_31, r_32, r_33], [V::from(0f32); BEAM_NODE_DOF / 2],

                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [r_11, r_12, r_13],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [r_21, r_22, r_23],
                [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
                [V::from(0f32); BEAM_NODE_DOF / 2], [r_31, r_32, r_33],
            ].concat(),
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


    fn dx_dr(x_1: V, x_2: V, r: V) -> V
    {
        BeamAuxFunctions::<V>::derivative_x(
            BeamAuxFunctions::<V>::power_func_x, x_1 * V::from(0.5f32), V::from(0f32),
            0) -
        BeamAuxFunctions::<V>::derivative_x(
            BeamAuxFunctions::<V>::power_func_x, x_1 * V::from(0.5f32), r, 1) +
        BeamAuxFunctions::<V>::derivative_x(
            BeamAuxFunctions::<V>::power_func_x, x_2 * V::from(0.5f32), V::from(0f32),
            0) +
        BeamAuxFunctions::<V>::derivative_x(
            BeamAuxFunctions::<V>::power_func_x, x_2 * V::from(0.5f32), r, 1)
    }


    fn jacobian(node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, FENode<V>>) -> V
    {
        let length = BeamAuxFunctions::length(node_1_number, node_2_number, nodes);
        let x_1 = V::from(-1f32) * length / V::from(2f32);
        let x_2 = length / V::from(2f32);
        BeamAuxFunctions::<V>::dx_dr(x_1, x_2, r)
    }


    fn inverse_jacobian(node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, FENode<V>>)
        -> V
    {
        V::from(1f32) / BeamAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    }


    pub fn determinant_of_jacobian(
        node_1_number: u32, 
        node_2_number: u32, r: V,
        nodes: &HashMap<u32, FENode<V>>
    ) 
        -> V
    {
        BeamAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    }


    pub fn h1_r(r: V) -> V
    {
        V::from(0.5f32) * (V::from(1f32) - r)
    }


    pub fn h2_r(r: V) -> V
    {
        V::from(0.5f32) * (V::from(1f32) + r)
    }


    fn dh1_dr(r: V) -> V
    {
        BeamAuxFunctions::<V>::derivative_x(
            BeamAuxFunctions::<V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) -
        BeamAuxFunctions::<V>::derivative_x(
            BeamAuxFunctions::<V>::power_func_x, V::from(0.5f32), r, 1)
    }


    fn dh2_dr(r: V) -> V
    {
        BeamAuxFunctions::<V>::derivative_x(
            BeamAuxFunctions::<V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) +
        BeamAuxFunctions::<V>::derivative_x(
            BeamAuxFunctions::<V>::power_func_x, V::from(0.5f32), r, 1)
    }


    pub fn strain_displacement_matrix_u(
        node_1_number: u32, 
        node_2_number: u32, 
        r: V, 
        nodes: &HashMap<u32, FENode<V>>
    ) 
        -> Result<Matrix<V>, String>
    {
        let elements = [
            BeamAuxFunctions::<V>::dh1_dr(r), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32),
            BeamAuxFunctions::<V>::dh2_dr(r), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32),
        ];
        let mut matrix = Matrix::create(
            1,
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &elements,
        );
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number, r, nodes);
        matrix.multiply_by_scalar(inverse_jacobian);
        Ok(matrix)
    }


    pub fn strain_displacement_matrix_v(
        node_1_number: u32,
        node_2_number: u32, 
        r: V,
        nodes: &HashMap<u32, FENode<V>>,
    ) 
        -> Result<Matrix<V>, String>
    {
        let lhs_elements = [
            V::from(0f32), BeamAuxFunctions::<V>::dh1_dr(r), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), BeamAuxFunctions::<V>::dh2_dr(r), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32),
        ];
        let mut lhs_matrix = Matrix::create(
            1,
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &lhs_elements,
        );
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number, r, nodes);
        lhs_matrix.multiply_by_scalar(inverse_jacobian);

        let rhs_elements = [
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), BeamAuxFunctions::<V>::h1_r(r),
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), BeamAuxFunctions::<V>::h2_r(r),
        ];

        let rhs_matrix = Matrix::create(
            1,
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &rhs_elements,
        );

        let matrix = lhs_matrix.subtract(&rhs_matrix)?;

        Ok(matrix)
    }


    pub fn strain_displacement_matrix_w(
        node_1_number: u32, 
        node_2_number: u32, 
        r: V,
        nodes: &HashMap<u32, FENode<V>>
    ) 
        -> Result<Matrix<V>, String>
    {
        let lhs_elements = [
            V::from(0f32), V::from(0f32), BeamAuxFunctions::<V>::dh1_dr(r),
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), BeamAuxFunctions::<V>::dh2_dr(r),
            V::from(0f32), V::from(0f32), V::from(0f32),
        ];
        let mut lhs_matrix = Matrix::create(
            1,
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &lhs_elements,
        );
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number, r, nodes);
        lhs_matrix.multiply_by_scalar(inverse_jacobian);

        let rhs_elements = [
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), BeamAuxFunctions::<V>::h1_r(r), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), BeamAuxFunctions::<V>::h2_r(r), V::from(0f32),
        ];
        let rhs_matrix = Matrix::create(
            1,
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &rhs_elements,
        );

        let matrix = lhs_matrix.subtract(&rhs_matrix)?;

        Ok(matrix)
    }


    pub fn strain_displacement_matrix_thu(
        node_1_number: u32, 
        node_2_number: u32, 
        r: V,
        nodes: &HashMap<u32, FENode<V>>,
    ) 
        -> Result<Matrix<V>, String>
    {
        let elements = [
            V::from(0f32), V::from(0f32), V::from(0f32),
            BeamAuxFunctions::<V>::dh1_dr(r), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32),
            BeamAuxFunctions::<V>::dh2_dr(r), V::from(0f32), V::from(0f32),
        ];
        let mut matrix = Matrix::create(
            1,
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &elements,
        );
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number, r, nodes);
        matrix.multiply_by_scalar(inverse_jacobian);
        Ok(matrix)
    }


    pub fn strain_displacement_matrix_thv(
        node_1_number: u32, 
        node_2_number: u32, 
        r: V, 
        nodes: &HashMap<u32, FENode<V>>
    ) 
        -> Result<Matrix<V>, String>
    {
        let elements = [
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), BeamAuxFunctions::<V>::dh1_dr(r), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), BeamAuxFunctions::<V>::dh2_dr(r), V::from(0f32),
        ];
        let mut matrix = Matrix::create(
            1,
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &elements,
        );
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number, r, nodes);
        matrix.multiply_by_scalar(inverse_jacobian);
        Ok(matrix)
    }


    pub fn strain_displacement_matrix_thw(
        node_1_number: u32, 
        node_2_number: u32, 
        r: V,
        nodes: &HashMap<u32, FENode<V>>,
    ) 
        -> Result<Matrix<V>, String>
    {
        let elements = [
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), BeamAuxFunctions::<V>::dh1_dr(r),
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), BeamAuxFunctions::<V>::dh2_dr(r),
        ];
        let mut matrix = Matrix::create(
            1,
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &elements,
        );
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number, r, nodes);
        matrix.multiply_by_scalar(inverse_jacobian);
        Ok(matrix)
    }


    pub fn local_stiffness_matrix(
        node_1_number: u32, 
        node_2_number: u32, 
        young_modulus: V,
        poisson_ratio: V, 
        area: V, 
        shear_factor: V, 
        i11: V, 
        i22: V,
        it: V, 
        alpha: V, 
        r: V, 
        nodes: &HashMap<u32, FENode<V>>,
    )
        -> Result<Matrix<V>, String>
    {
        let mut lhs_matrix_u = BeamAuxFunctions::strain_displacement_matrix_u(
                node_1_number, node_2_number, r, nodes,
            )?
            .transpose();

        let rhs_matrix_u = BeamAuxFunctions::strain_displacement_matrix_u(
            node_1_number, node_2_number, r, nodes,
        )?;

        let mut matrix_u = lhs_matrix_u
            .multiply(&rhs_matrix_u)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {e}"))?
            .multiply_by_scalar(
                young_modulus * 
                area * 
                BeamAuxFunctions::determinant_of_jacobian(node_1_number, node_2_number, r, nodes) * 
                alpha
            );

        let mut lhs_matrix_v = BeamAuxFunctions::strain_displacement_matrix_v(
                node_1_number, node_2_number, r, nodes,
            )
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .transpose();

        let rhs_matrix_v = BeamAuxFunctions::strain_displacement_matrix_v(
                node_1_number, node_2_number, r, nodes,
            )
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;

        let shear_modulus = young_modulus / (V::from(2f32) * (V::from(1f32) + poisson_ratio));
        let mut matrix_v = lhs_matrix_v
            .multiply(&rhs_matrix_v)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .multiply_by_scalar(
                shear_modulus * 
                area *
                shear_factor * 
                BeamAuxFunctions::determinant_of_jacobian(node_1_number, node_2_number, r, nodes) * 
                alpha
            );

        let mut lhs_matrix_w = BeamAuxFunctions::strain_displacement_matrix_w(
                node_1_number, node_2_number, r, nodes,
            )
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .transpose();

        let rhs_matrix_w = BeamAuxFunctions::strain_displacement_matrix_w(
                node_1_number, node_2_number, r, nodes,
            )
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;

        let mut matrix_w = lhs_matrix_w.multiply(&rhs_matrix_w)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .multiply_by_scalar(
                shear_modulus * 
                area * 
                shear_factor *
                BeamAuxFunctions::determinant_of_jacobian(node_1_number, node_2_number, r, nodes) *
                alpha
            );

        let mut lhs_matrix_thu = BeamAuxFunctions::strain_displacement_matrix_thu(
                node_1_number, node_2_number, r, nodes
            )?
            .transpose();

        let rhs_matrix_thu = BeamAuxFunctions::strain_displacement_matrix_thu(
            node_1_number, node_2_number, r, nodes,
        )?;

        let mut matrix_thu = lhs_matrix_thu
            .multiply(&rhs_matrix_thu)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .multiply_by_scalar(
                shear_modulus * 
                it *
                BeamAuxFunctions::determinant_of_jacobian(node_1_number, node_2_number, r, nodes) * 
                alpha
            );

        let mut lhs_matrix_thv = BeamAuxFunctions::strain_displacement_matrix_thv(
                node_1_number, node_2_number, r, nodes,
            )?
            .transpose();

        let rhs_matrix_thv = BeamAuxFunctions::strain_displacement_matrix_thv(
            node_1_number, node_2_number, r, nodes,
        )?;

        let mut matrix_thv = lhs_matrix_thv
            .multiply(&rhs_matrix_thv)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .multiply_by_scalar(
                young_modulus *
                i22 *
                BeamAuxFunctions::determinant_of_jacobian(node_1_number, node_2_number, r, nodes) * 
                alpha
            );

        let mut lhs_matrix_thw = BeamAuxFunctions::strain_displacement_matrix_thw(
                node_1_number, node_2_number, r, nodes
            )?
            .transpose();

        let rhs_matrix_thw = BeamAuxFunctions::strain_displacement_matrix_thw(
            node_1_number, node_2_number, r, nodes,
        )?;

        let mut matrix_thw = lhs_matrix_thw
            .multiply(&rhs_matrix_thw)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .multiply_by_scalar(
                young_modulus * 
                i11 *
                BeamAuxFunctions::determinant_of_jacobian(node_1_number, node_2_number, r, nodes) * 
                alpha
            );

        let matrix = matrix_u
            .add(&matrix_v)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .add(&matrix_w)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .add(&matrix_thu)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .add(&matrix_thv)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?
            .add(&matrix_thw)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;

        Ok(matrix)
    }


    pub fn compose_node_dof_parameters(node_number: u32) -> Result<Vec<DOFParameterData>, String>
    {
        let mut node_dof_parameters = Vec::new();
        for dof in 0..BEAM_NODE_DOF
        {
            let dof_parameter =
                GlobalDOFParameter::iterator().nth(dof)
                    .ok_or("Beam2n2ipT: Could not find dof parameter!")?;
            let dof_parameter = DOFParameterData::create(node_number,
                *dof_parameter);
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
                let value = column_matrix.get_element_value(&Position(row, column))?;
                values.push(*value);
                column += 1;
            }
            row += 1;
        }
        Ok(values)
    }


    pub fn find_principal_moments_of_inertia(i11_init: V, i22_init: V, i12_init: V) -> (V, V, V)
    {
        // if compare_with_tolerance(i12_init) == V::from(0f32)
        // {
        //     return 
        // }
        // let mut angle = if i11_init == i22_init
        //     {
        //         V::from(0f32)
        //     }
        //     else
        //     {
        //         (V::from(2f32) * i12_init / (i22_init - i11_init)).my_atan() /
        //         V::from(2f32)
        //     };
        let mut angle = (V::from(2f32) * i12_init / (i22_init - i11_init)).my_atan() / V::from(2f32);

        let mut i11 = i11_init * (angle.my_cos()).my_powi(2) +
            i22_init * (angle.my_sin()).my_powi(2) -
            i12_init * (V::from(2f32) * angle).my_sin();

        let mut i22 = i11_init * (angle.my_sin()).my_powi(2) +
            i22_init * (angle.my_cos()).my_powi(2) +
            i12_init * (V::from(2f32) * angle).my_sin();

        let mut i = 1;
        while i11 < i22
        {
            angle = ((V::from(2f32) * i12_init / (i22_init - i11_init)).my_atan() +
                V::from(PI) * V::from(i as f32)) / V::from(2f32);
            i11 = i11_init * (angle.my_cos()).my_powi(2) +
                i22_init * (angle.my_sin()).my_powi(2) -
                i12_init * (V::from(2f32) * angle).my_sin();
            i22 = i11_init * (angle.my_sin()).my_powi(2) +
                i22_init * (angle.my_cos()).my_powi(2) +
                i12_init * (V::from(2f32) * angle).my_sin();
            i += 1;
        }

        (i11, i22, angle)
    }
}
