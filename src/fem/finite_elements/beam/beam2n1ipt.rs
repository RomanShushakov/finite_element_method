use std::hash::Hash;
use std::fmt::Debug;
use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};
use std::collections::HashMap;
use std::f32::consts::PI;

use extended_matrix::basic_matrix::basic_matrix::MatrixElementPosition;
use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::extract_element_value;

use crate::fem::finite_elements::finite_element::FiniteElementTrait;
use crate::fem::finite_elements::fe_node::FENode;

use crate::fem::global_analysis::fe_stiffness::{StiffnessGroup, StiffnessType};
use crate::fem::global_analysis::fe_dof_parameter_data::{GlobalDOFParameter, DOFParameterData};
use crate::fem::global_analysis::fe_global_analysis_result::Displacements;

use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::element_analysis::fe_element_analysis_result::{ElementAnalysisData, ElementForces};
use crate::fem::element_analysis::fe_stress_strain_components::StressStrainComponent;

use crate::fem::finite_elements::functions::compare_with_tolerance;

use crate::my_float::MyFloatTrait;

use crate::fem::finite_elements::beam::consts::{ BEAM_NODE_DOF, BEAM2N1IPT_NODES_NUMBER };


struct BeamAuxFunctions<T, V>(T, V);


impl<T, V> BeamAuxFunctions<T, V>
    where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
             From<u8> + 'static,
          V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + From<f32> +
             Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
             MyFloatTrait + PartialOrd + MyFloatTrait<Other = V> + 'static,
{
    fn length(node_1_number: T, node_2_number: T, nodes: &HashMap<T, FENode<V>>) -> V
    {
        let node_1 = nodes.get(&node_1_number).unwrap();
        let node_2 = nodes.get(&node_2_number).unwrap();

        ((node_1.x() - node_2.x()).my_powi(2) + (node_1.y() - node_2.y()).my_powi(2) +
        (node_1.z() - node_2.z()).my_powi(2)).my_sqrt()
    }


    fn nodes_number() -> T
    {
        let mut n = T::from(0u8);
        (0..BEAM2N1IPT_NODES_NUMBER).for_each(|_| n += T::from(1u8));
        n
    }


    fn node_dof() -> T
    {
        let mut n = T::from(0u8);
        (0..BEAM_NODE_DOF).for_each(|_| n += T::from(1u8));
        n
    }


    fn find_components_of_line_a_perpendicular_to_line_b(line_a: &[V; 3],
        line_b: &[V; 3], tolerance: V) -> Result<[V; 3], String>
    {
        let a_x = V::from(-1f32) * line_a[0];
        let a_y = V::from(-1f32) * line_a[1];
        let a_z = V::from(-1f32) * line_a[2];

        let a = ExtendedMatrix::create(T::from(3u8),
            T::from(1u8), vec![a_x, a_y, a_z], tolerance);

        let b_x = line_b[0];
        let b_y = line_b[1];
        let b_z = line_b[2];

        let norm = V::from(1f32) / (b_x.my_powi(2) + b_y.my_powi(2) + b_z.my_powi(2));

        let mut coeff_matrix = ExtendedMatrix::create(T::from(3u8),
            T::from(3u8), vec![
                V::from(-1f32) * b_z * b_z - b_y * b_y, b_x * b_y, b_x * b_z,
                b_y * b_x, V::from(-1f32) * b_x * b_x - b_z * b_z,	b_y * b_z,
                b_z * b_x,	b_z * b_y, V::from(-1f32) * b_y * b_y - b_x * b_x,
            ], tolerance);

        coeff_matrix.multiply_by_number(norm);

        let components_of_line_a_perpendicular_to_line_b_matrix = coeff_matrix
            .multiply_by_matrix(&a)?;

        let components_of_line_a_perpendicular_to_line_b_all_values =
            components_of_line_a_perpendicular_to_line_b_matrix.extract_all_elements_values();

        let a_perpendicular_to_b_x = extract_element_value(T::from(0u8), T::from(0u8),
            &components_of_line_a_perpendicular_to_line_b_all_values);

        let a_perpendicular_to_b_y = extract_element_value(T::from(1u8), T::from(0u8),
            &components_of_line_a_perpendicular_to_line_b_all_values);

        let a_perpendicular_to_b_z = extract_element_value(T::from(2u8), T::from(0u8),
            &components_of_line_a_perpendicular_to_line_b_all_values);

        Ok([a_perpendicular_to_b_x, a_perpendicular_to_b_y, a_perpendicular_to_b_z])
    }


    fn rotation_matrix(node_1_number: T, node_2_number: T, local_axis_1_direction: &[V; 3],
        angle: V, tolerance: V, nodes: &HashMap<T, FENode<V>>)
        -> Result<ExtendedMatrix<T, V>, String>
    {
        let node_1 = nodes.get(&node_1_number).unwrap();
        let node_2 = nodes.get(&node_2_number).unwrap();

        let x = node_2.x() - node_1.x();
        let y = node_2.y() - node_1.y();
        let z = node_2.z() - node_1.z();

        let length = BeamAuxFunctions::<T, V>::length(node_1_number, node_2_number, nodes);

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
            rotation_axis_coord_z = V::from(-1f32) * y * length;
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

        let components_projection_of_beam_section_orientation_vector =
            BeamAuxFunctions::<T, V>::find_components_of_line_a_perpendicular_to_line_b(
                local_axis_1_direction, &[x, y, z],
                tolerance
            )?;

        let interim_rotation_matrix = ExtendedMatrix::create(
            3, 3, vec![q_11, q_12, q_13, q_21, q_22, q_23, q_31,
            q_32, q_33], tolerance);

        let projection_of_beam_section_orientation = ExtendedMatrix::create(
            3, 1,
            components_projection_of_beam_section_orientation_vector.to_vec(),
            tolerance);

        let transformed_projection_of_beam_section_orientation =
            interim_rotation_matrix.multiply_by_matrix(&projection_of_beam_section_orientation)?;

        let all_values_of_transformed_projection_of_beam_section_orientation =
            transformed_projection_of_beam_section_orientation.extract_all_elements_values();

        let transformed_projection_of_beam_section_orientation_x = extract_element_value(0,
            0, &all_values_of_transformed_projection_of_beam_section_orientation);

        let transformed_projection_of_beam_section_orientation_y = extract_element_value(1,
            0, &all_values_of_transformed_projection_of_beam_section_orientation);

        let transformed_projection_of_beam_section_orientation_z = extract_element_value(2,
            0, &all_values_of_transformed_projection_of_beam_section_orientation);

        let angle_between_beam_section_local_axis_1_direction_and_horizon =
            (V::from(-1f32) * transformed_projection_of_beam_section_orientation_z /
            (transformed_projection_of_beam_section_orientation_x.my_powi(2) +
            transformed_projection_of_beam_section_orientation_y.my_powi(2) +
            transformed_projection_of_beam_section_orientation_z.my_powi(2))
                .my_sqrt()
            ).my_acos();

        let total_angle = angle + angle_between_beam_section_local_axis_1_direction_and_horizon;

        let c_x = compare_with_tolerance(x /length, tolerance);
        let c_y = compare_with_tolerance(y / length, tolerance);
        let c_z = compare_with_tolerance(z / length, tolerance);
        let c_xz = compare_with_tolerance(
            (c_x.my_powi(2) + c_z.my_powi(2)).my_sqrt(), tolerance);

        let c = total_angle.my_cos();
        let s = total_angle.my_sin();

        let r_11 = if c_xz != V::from(0f32) { c_x } else { V::from(0f32) };
        let r_12 = c_y;
        let r_13 = if c_xz != V::from(0f32) { c_z } else { V::from(0f32) };
        let r_21 = if c_xz != V::from(0f32)
            { (V::from(-1f32) * c_x * c_y * c - c_z * s) / c_xz } else { V::from(-1f32) * c_y * c };
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
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            vec![
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
            tolerance,
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
        BeamAuxFunctions::<T, V>::derivative_x(
            BeamAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.5f32), V::from(0f32),
            0) -
        BeamAuxFunctions::<T, V>::derivative_x(
            BeamAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.5f32), r, 1) +
        BeamAuxFunctions::<T, V>::derivative_x(
            BeamAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.5f32), V::from(0f32),
            0) +
        BeamAuxFunctions::<T, V>::derivative_x(
            BeamAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.5f32), r, 1)
    }


    fn jacobian(node_1_number: T, node_2_number: T, r: V, nodes: &HashMap<T, FENode<V>>) -> V
    {
        let length = BeamAuxFunctions::length(node_1_number, node_2_number, nodes);
        let x_1 = V::from(-1f32) * length / V::from(2f32);
        let x_2 = length / V::from(2f32);
        BeamAuxFunctions::<T, V>::dx_dr(x_1, x_2, r)
    }


    fn inverse_jacobian(node_1_number: T, node_2_number: T, r: V, nodes: &HashMap<T, FENode<V>>) -> V
    {
        V::from(1f32) / BeamAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    }


    fn determinant_of_jacobian(node_1_number: T, node_2_number: T, r: V,
        nodes: &HashMap<T, FENode<V>>) -> V
    {
        BeamAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    }


    fn h1_r(r: V) -> V
    {
        V::from(0.5f32) * (V::from(1f32) - r)
    }


    fn h2_r(r: V) -> V
    {
        V::from(0.5f32) * (V::from(1f32) + r)
    }


    fn dh1_dr(r: V) -> V
    {
        BeamAuxFunctions::<T, V>::derivative_x(
            BeamAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) -
        BeamAuxFunctions::<T, V>::derivative_x(
            BeamAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), r, 1)
    }


    fn dh2_dr(r: V) -> V
    {
        BeamAuxFunctions::<T, V>::derivative_x(
            BeamAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) +
        BeamAuxFunctions::<T, V>::derivative_x(
            BeamAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), r, 1)
    }


    fn strain_displacement_matrix_u(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> ExtendedMatrix<T, V>
    {
        let elements = vec![
                BeamAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
                BeamAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
            ];
        let mut matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            elements, tolerance);
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        matrix.multiply_by_number(inverse_jacobian);
        matrix
    }


    fn strain_displacement_matrix_v(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let lhs_elements = vec![
                V::from(0f32), BeamAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), BeamAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
            ];
        let mut lhs_matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            lhs_elements, tolerance);
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        lhs_matrix.multiply_by_number(inverse_jacobian);

        let rhs_elements = vec![
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), BeamAuxFunctions::<T, V>::h1_r(r),
            V::from(0f32), V::from(0f32), V::from(0f32),
            V::from(0f32), V::from(0f32), BeamAuxFunctions::<T, V>::h2_r(r),
        ];

        let rhs_matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            rhs_elements, tolerance);

        let matrix = lhs_matrix.subtract_matrix(&rhs_matrix)?;

        Ok(matrix)
    }


    fn strain_displacement_matrix_w(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let lhs_elements = vec![
                V::from(0f32), V::from(0f32), BeamAuxFunctions::<T, V>::dh1_dr(r),
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), BeamAuxFunctions::<T, V>::dh2_dr(r),
                V::from(0f32), V::from(0f32), V::from(0f32),
            ];
        let mut lhs_matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            lhs_elements, tolerance);
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        lhs_matrix.multiply_by_number(inverse_jacobian);

        let rhs_elements = vec![
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), BeamAuxFunctions::<T, V>::h1_r(r), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), BeamAuxFunctions::<T, V>::h2_r(r), V::from(0f32),
            ];
        let rhs_matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            rhs_elements, tolerance);

        let matrix = lhs_matrix.subtract_matrix(&rhs_matrix)?;

        Ok(matrix)
    }


    fn strain_displacement_matrix_thu(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> ExtendedMatrix<T, V>
    {
        let elements = vec![
                V::from(0f32), V::from(0f32), V::from(0f32),
                BeamAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
                BeamAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32), V::from(0f32),
            ];
        let mut matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            elements, tolerance);
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        matrix.multiply_by_number(inverse_jacobian);
        matrix
    }


    fn strain_displacement_matrix_thv(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> ExtendedMatrix<T, V>
    {
        let elements = vec![
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), BeamAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), BeamAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32),
            ];
        let mut matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            elements, tolerance);
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        matrix.multiply_by_number(inverse_jacobian);
        matrix
    }


    fn strain_displacement_matrix_thw(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> ExtendedMatrix<T, V>
    {
        let elements = vec![
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), BeamAuxFunctions::<T, V>::dh1_dr(r),
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), BeamAuxFunctions::<T, V>::dh2_dr(r),
            ];
        let mut matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            elements, tolerance);
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        matrix.multiply_by_number(inverse_jacobian);
        matrix
    }


    fn local_stiffness_matrix(node_1_number: T, node_2_number: T, young_modulus: V,
        poisson_ratio: V, area: V, shape_factor: V, it: V, i11: V, i22: V,
        alpha: V, r: V, tolerance: V, nodes: &HashMap<T, FENode<V>>)
        -> Result<ExtendedMatrix<T, V>, String>
    {
        // let f = |data: &str| println!("{}", data);

        let mut lhs_matrix_u = BeamAuxFunctions::strain_displacement_matrix_u(
            node_1_number, node_2_number, r, tolerance, nodes);
        lhs_matrix_u.transpose();
        let rhs_matrix_u = BeamAuxFunctions::strain_displacement_matrix_u(
            node_1_number, node_2_number, r, tolerance, nodes);
        let mut matrix_u = lhs_matrix_u.multiply_by_matrix(&rhs_matrix_u)
            .map_err(|e|
                format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let coeff_u = young_modulus * area;
        matrix_u.multiply_by_number(coeff_u);
        matrix_u.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

        // println!("Matrix u");
        // matrix_u.show_matrix(f);
        // println!();

        let mut lhs_matrix_v = BeamAuxFunctions::strain_displacement_matrix_v(
            node_1_number, node_2_number, r, tolerance, nodes)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        lhs_matrix_v.transpose();
        let rhs_matrix_v = BeamAuxFunctions::strain_displacement_matrix_v(
            node_1_number, node_2_number, r, tolerance, nodes)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let mut matrix_v = lhs_matrix_v.multiply_by_matrix(&rhs_matrix_v)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let shear_modulus = young_modulus / (V::from(2f32) * (V::from(1f32) + poisson_ratio));
        let coeff_v_and_w = shear_modulus * area * shape_factor;
        matrix_v.multiply_by_number(coeff_v_and_w);
        matrix_v.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

        // println!("Matrix v");
        // matrix_v.show_matrix(f);
        // println!();

        let mut lhs_matrix_w = BeamAuxFunctions::strain_displacement_matrix_w(
            node_1_number, node_2_number, r, tolerance, nodes)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        lhs_matrix_w.transpose();
        let rhs_matrix_w = BeamAuxFunctions::strain_displacement_matrix_w(
            node_1_number, node_2_number, r, tolerance, nodes)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let mut matrix_w = lhs_matrix_w.multiply_by_matrix(&rhs_matrix_w)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        matrix_w.multiply_by_number(coeff_v_and_w);
        matrix_w.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

        // println!("Matrix w");
        // matrix_w.show_matrix(f);
        // println!();

        let mut lhs_matrix_thu = BeamAuxFunctions::strain_displacement_matrix_thu(
            node_1_number, node_2_number, r, tolerance, nodes);
        lhs_matrix_thu.transpose();
        let rhs_matrix_thu = BeamAuxFunctions::strain_displacement_matrix_thu(
            node_1_number, node_2_number, r, tolerance, nodes);
        let mut matrix_thu = lhs_matrix_thu.multiply_by_matrix(
            &rhs_matrix_thu)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let coeff_thu = shear_modulus * it;
        matrix_thu.multiply_by_number(coeff_thu);
        matrix_thu.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

        // println!("Matrix thu");
        // matrix_thu.show_matrix(f);
        // println!();

        let mut lhs_matrix_thv = BeamAuxFunctions::strain_displacement_matrix_thv(
            node_1_number, node_2_number, r, tolerance, nodes);
        lhs_matrix_thv.transpose();
        let rhs_matrix_thv = BeamAuxFunctions::strain_displacement_matrix_thv(
            node_1_number, node_2_number, r, tolerance, nodes);
        let mut matrix_thv = lhs_matrix_thv.multiply_by_matrix(
            &rhs_matrix_thv)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let coeff_thv = young_modulus * i22;
        matrix_thv.multiply_by_number(coeff_thv);
        matrix_thv.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

        // println!("Matrix thv");
        // matrix_thv.show_matrix(f);
        // println!();

        let mut lhs_matrix_thw = BeamAuxFunctions::strain_displacement_matrix_thw(
            node_1_number, node_2_number, r, tolerance, nodes);
        lhs_matrix_thw.transpose();
        let rhs_matrix_thw = BeamAuxFunctions::strain_displacement_matrix_thw(
            node_1_number, node_2_number, r, tolerance, nodes);
        let mut matrix_thw = lhs_matrix_thw.multiply_by_matrix(
            &rhs_matrix_thw)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let coeff_thw = young_modulus * i11;
        matrix_thw.multiply_by_number(coeff_thw);
        matrix_thw.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

        // println!("Matrix thw");
        // matrix_thw.show_matrix(f);
        // println!();

        let matrix = ((((matrix_u.add_matrix(&matrix_v)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?).add_matrix(&matrix_w)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?).add_matrix(&matrix_thu)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?).add_matrix(&matrix_thv)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?).add_matrix(&matrix_thw)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;

        Ok(matrix)
    }


    fn compose_node_dof_parameters<'a>(node_number: T)
        -> Result<Vec<DOFParameterData<T>>, &'a str>
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


    fn extract_column_matrix_values(column_matrix: &ExtendedMatrix<T, V>) -> Vec<V>
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


struct IntegrationPoint<V>
{
    r: V,
    weight: V,
}


struct State<T, V>
{
    rotation_matrix: ExtendedMatrix<T, V>,
    integration_points: Vec<IntegrationPoint<V>>,
    local_stiffness_matrix: ExtendedMatrix<T, V>,
    nodes_dof_parameters_global: Vec<DOFParameterData<T>>,
}


impl<T, V> State<T, V>
{
    fn create(rotation_matrix: ExtendedMatrix<T, V>, integration_points: Vec<IntegrationPoint<V>>,
        local_stiffness_matrix: ExtendedMatrix<T, V>,
        nodes_dof_parameters_global: Vec<DOFParameterData<T>>) -> Self
    {
        State { rotation_matrix, integration_points, local_stiffness_matrix,
        nodes_dof_parameters_global }
    }
}


pub struct Beam2n1ipT<T, V>
{
    node_1_number: T,
    node_2_number: T,
    young_modulus: V,
    poisson_ratio: V,
    area: V,
    i11: V,
    i22: V,
    angle: V,
    shape_factor: V,
    it: V,
    local_axis_1_direction: [V; 3],
    state: State<T, V>,
}


impl<T, V> Beam2n1ipT<T, V>
    where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
             From<u8> + 'static,
          V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + From<f32> + Add<Output = V> +
             Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
             MyFloatTrait + PartialOrd + MyFloatTrait<Other = V> + 'static
{
    pub fn create(node_1_number: T, node_2_number: T, young_modulus: V, poisson_ratio: V, area: V,
        i11_init: V, i22_init: V, i12_init: V, it: V, shape_factor: V,
        local_axis_1_direction: [V; 3], tolerance: V, nodes: &HashMap<T, FENode<V>>)
        -> Result<Self, String>
    {
        let mut angle = if i11_init == i22_init
            {
                V::from(0f32)
            }
            else
            {
                (V::from(2f32) * i12_init / (i22_init - i11_init)).my_atan() /
                V::from(2f32)
            };

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

        let integration_point_1 = IntegrationPoint {
            r: V::from(0f32), weight: V::from(2f32) };

        let rotation_matrix =
            BeamAuxFunctions::<T, V>::rotation_matrix(node_1_number, node_2_number,
                &local_axis_1_direction, angle, tolerance, nodes)?;

        let integration_points = vec![integration_point_1];

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (BEAM2N1IPT_NODES_NUMBER * BEAM_NODE_DOF).pow(2)], tolerance);

        for integration_point in &integration_points
        {
            let matrix = BeamAuxFunctions::local_stiffness_matrix(
                node_1_number, node_2_number, young_modulus, poisson_ratio, area, shape_factor,
                it, i11, i22, integration_point.weight, integration_point.r, tolerance,
                nodes)?;

            local_stiffness_matrix = local_stiffness_matrix.add_matrix(&matrix)
                .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        let mut nodes_dof_parameters =
            BeamAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;

        let node_2_dof_parameters =
            BeamAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;

        nodes_dof_parameters.extend(node_2_dof_parameters);

        let state = State::create(rotation_matrix, integration_points,
            local_stiffness_matrix, nodes_dof_parameters);

        Ok(Beam2n1ipT { node_1_number, node_2_number, young_modulus, poisson_ratio,
            area, i11, i22, angle, shape_factor, it, local_axis_1_direction, state })
    }


    fn extract_local_displacements(&self, global_displacements: &Displacements<T, V>, tolerance: V)
        -> Result<ExtendedMatrix<T, V>, String>
    {
        let mut element_global_displacements_values = Vec::new();
        for lhs_dof_parameter_data in &self.state.nodes_dof_parameters_global
        {
            if let Some(position) = global_displacements.dof_parameters_data()
                .iter()
                .position(|rhs_dof_parameter_data|
                    rhs_dof_parameter_data == lhs_dof_parameter_data)
            {
                let displacement_value = global_displacements.displacements_values()[position];
                element_global_displacements_values.push(displacement_value);
            }
            else
            {
                element_global_displacements_values.push(V::from(0f32));
            }
        }

        let mut m = 0usize;
        let mut rows_number = T::from(0u8);
        while m < self.state.nodes_dof_parameters_global.len()
        {
            m += 1usize;
            rows_number += T::from(1u8);
        }

        let element_global_displacements = ExtendedMatrix::create(rows_number,
            T::from(1u8), element_global_displacements_values,
            tolerance);

        let element_local_displacements =
            self.state.rotation_matrix.multiply_by_matrix(&element_global_displacements)?;
        Ok(element_local_displacements)
    }
}


impl<T, V> FiniteElementTrait<T, V> for Beam2n1ipT<T, V>
    where T: Copy + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> +
             Mul<Output = T> + Eq + Hash + Debug + SubAssign + PartialOrd + AddAssign +
             From<u8> + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             Into<f64> + SubAssign + AddAssign + MulAssign + PartialEq + Debug +
             MyFloatTrait + PartialOrd + From<f32> + MyFloatTrait<Other = V> + 'static,
{
    fn update(&mut self, nodes_numbers: Vec<T>, properties: Vec<V>, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        let node_1_number = nodes_numbers[0];

        let node_2_number = nodes_numbers[1];

        let young_modulus = properties[0];

        let poisson_ratio = properties[1];

        let area = properties[2];

        let i11_init = properties[3];

        let i22_init = properties[4];

        let i12_init = properties[5];

        let mut angle = if i11_init == i22_init
            {
                V::from(0f32)
            }
            else
            {
                (V::from(2f32) * i12_init / (i22_init - i11_init)).my_atan() /
                V::from(2f32)
            };

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

        let it = properties[6];

        let shape_factor = properties[7];

        let local_axis_1_direction = [properties[8], properties[9], properties[10]];

        let rotation_matrix =
            BeamAuxFunctions::rotation_matrix(node_1_number, node_2_number,
                &local_axis_1_direction, angle, tolerance, nodes)?;

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (BEAM2N1IPT_NODES_NUMBER * BEAM_NODE_DOF).pow(2)], tolerance);

        for integration_point in &self.state.integration_points
        {
            let matrix = BeamAuxFunctions::local_stiffness_matrix(
                node_1_number, node_2_number, young_modulus, poisson_ratio, area, shape_factor,
                it, i11, i22, integration_point.weight, integration_point.r, tolerance,
                nodes)?;
            local_stiffness_matrix = local_stiffness_matrix.add_matrix(&matrix)
                .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        let mut nodes_dof_parameters =
            BeamAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;

        let node_2_dof_parameters =
            BeamAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;

        nodes_dof_parameters.extend(node_2_dof_parameters);

        self.node_1_number = node_1_number;
        self.node_2_number = node_2_number;
        self.young_modulus = young_modulus;
        self.poisson_ratio = poisson_ratio;
        self.area = area;
        self.i11 = i11;
        self.i22 = i22;
        self.angle = angle;
        self.shape_factor = shape_factor;
        self.it = it;
        self.local_axis_1_direction = local_axis_1_direction;
        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        self.state.nodes_dof_parameters_global = nodes_dof_parameters;

        Ok(())
    }


    fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
    {
        let mut interim_matrix = self.state.rotation_matrix.clone();
        interim_matrix.transpose();
        if let Ok(matrix) =
        interim_matrix.multiply_by_matrix(&self.state.local_stiffness_matrix)
        {
            if let Ok(matrix) =
            matrix.multiply_by_matrix(&self.state.rotation_matrix)
            {
                return Ok(matrix);
            }
        }
        Err("Beam2n2ipT: Stiffness matrix cannot be extracted!")
    }


    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>
    {
        let (rows_number, columns_number) =
            (BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
             BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof());

        let mut positions_kuu_1_1 = Vec::new();
        let mut positions_kuth_1_1 = Vec::new();
        let mut positions_kthu_1_1 = Vec::new();
        let mut positions_kthth_1_1 = Vec::new();
        let mut positions_kuu_1_2 = Vec::new();
        let mut positions_kuth_1_2 = Vec::new();
        let mut positions_kthu_1_2 = Vec::new();
        let mut positions_kthth_1_2 = Vec::new();
        let mut positions_kuu_2_1 = Vec::new();
        let mut positions_kuth_2_1 = Vec::new();
        let mut positions_kthu_2_1 = Vec::new();
        let mut positions_kthth_2_1 = Vec::new();
        let mut positions_kuu_2_2 = Vec::new();
        let mut positions_kuth_2_2 = Vec::new();
        let mut positions_kthu_2_2 = Vec::new();
        let mut positions_kthth_2_2 = Vec::new();

        let mut i = T::from(0u8);
        while i < rows_number * columns_number
        {
            let position = MatrixElementPosition::create(
                i / columns_number, i % columns_number);

            let row = i / columns_number;
            let column = i % columns_number;

            if row < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)
            {
                positions_kuu_1_1.push(position);
            }
            else if row < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof()
            {
                positions_kuth_1_1.push(position);
            }
            else if row < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() &&
                column < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kuu_1_2.push(position);
            }
            else if row < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kuth_1_2.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                row < BeamAuxFunctions::<T, V>::node_dof() &&
                column < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)
            {
                positions_kthu_1_1.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                row < BeamAuxFunctions::<T, V>::node_dof() &&
                column >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof()
            {
                positions_kthth_1_1.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                row < BeamAuxFunctions::<T, V>::node_dof() &&
                column >= BeamAuxFunctions::<T, V>::node_dof() &&
                column < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kthu_1_2.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                row < BeamAuxFunctions::<T, V>::node_dof() &&
                column >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kthth_1_2.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() &&
                row < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)
            {
                positions_kuu_2_1.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() &&
                row < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof()
            {
                positions_kuth_2_1.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() &&
                row < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() &&
                column < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kuu_2_2.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() &&
                row < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kuth_2_2.push(position);
            }
            else if row >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)
            {
                positions_kthu_2_1.push(position);
            }
            else if row >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof()
            {
                positions_kthth_2_1.push(position);
            }
            else if row >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() &&
                column < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kthu_2_2.push(position);
            }
            else
            {
                positions_kthth_2_2.push(position);
            }
            i += T::from(1u8);
        }

        vec![StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                number_2: self.node_1_number, positions: positions_kuu_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_1_number,
                number_2: self.node_1_number, positions: positions_kuth_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                number_2: self.node_2_number, positions: positions_kuu_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_1_number,
                number_2: self.node_2_number, positions: positions_kuth_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_1_number,
                number_2: self.node_1_number, positions: positions_kthu_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_1_number,
                number_2: self.node_1_number, positions: positions_kthth_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_1_number,
                number_2: self.node_2_number, positions: positions_kthu_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_1_number,
                number_2: self.node_2_number, positions: positions_kthth_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                number_2: self.node_1_number, positions: positions_kuu_2_1 },
             StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_2_number,
                number_2: self.node_1_number, positions: positions_kuth_2_1 },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                number_2: self.node_2_number, positions: positions_kuu_2_2 },
             StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_2_number,
                number_2: self.node_2_number, positions: positions_kuth_2_2 },
             StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_2_number,
                number_2: self.node_1_number, positions: positions_kthu_2_1 },
             StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_2_number,
                number_2: self.node_1_number, positions: positions_kthth_2_1 },
             StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_2_number,
                number_2: self.node_2_number, positions: positions_kthu_2_2 },
             StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_2_number,
                number_2: self.node_2_number, positions: positions_kthth_2_2 },
        ]
    }


    fn node_belong_element(&self, node_number: T) -> bool
    {
        self.node_1_number == node_number || self.node_2_number == node_number
    }


    fn refresh(&mut self, tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        let rotation_matrix =
            BeamAuxFunctions::rotation_matrix(self.node_1_number, self.node_2_number,
                &self.local_axis_1_direction, self.angle, tolerance, nodes)?;

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (BEAM2N1IPT_NODES_NUMBER * BEAM_NODE_DOF).pow(2)], tolerance);

        for integration_point in self.state.integration_points.iter()
        {
            let matrix = BeamAuxFunctions::local_stiffness_matrix(
                self.node_1_number, self.node_2_number, self.young_modulus, self.poisson_ratio,
                self.area, self.shape_factor, self.it, self.i11, self.i22,
                integration_point.weight, integration_point.r, tolerance, nodes)?;
            local_stiffness_matrix = local_stiffness_matrix.add_matrix(&matrix)
                .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        Ok(())
    }


    fn nodes_numbers_same(&self, nodes_numbers: Vec<T>) -> bool
    {
        (nodes_numbers[0] == self.node_1_number && nodes_numbers[1] == self.node_2_number) ||
        (nodes_numbers[0] == self.node_2_number && nodes_numbers[1] == self.node_1_number)
    }


    fn extract_element_analysis_data(&self, global_displacements: &Displacements<T, V>,
        tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<ElementAnalysisData<V>, String>
    {
        let element_local_displacements =
            self.extract_local_displacements(global_displacements, tolerance)?;

        let r = self.state.integration_points[0].r;
        let mut forces_values = Vec::new();
        let mut forces_components = Vec::new();

        let strain_displacement_matrix_u =
            BeamAuxFunctions::strain_displacement_matrix_u(
                self.node_1_number, self.node_2_number, r, tolerance, nodes);
        let strains_matrix_u =
            strain_displacement_matrix_u.multiply_by_matrix(&element_local_displacements)?;
        let force_x = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_u)[0] * self.area * self.young_modulus;
        forces_components.push(ForceComponent::ForceX);
        forces_values.push(force_x);

        let strain_displacement_matrix_v =
            BeamAuxFunctions::strain_displacement_matrix_v(
                self.node_1_number, self.node_2_number, r, tolerance, nodes)?;
        let strains_matrix_v =
            strain_displacement_matrix_v.multiply_by_matrix(&element_local_displacements)?;
        let shear_modulus = self.young_modulus /
            (V::from(2f32) * (V::from(1f32) + self.poisson_ratio));
        let shear_y = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_v)[0] * shear_modulus * self.area * self.shape_factor;
        forces_components.push(ForceComponent::ForceY);
        forces_values.push(shear_y);

        let strain_displacement_matrix_w =
            BeamAuxFunctions::strain_displacement_matrix_w(
                self.node_1_number, self.node_2_number, r, tolerance, nodes)?;
        let strains_matrix_w =
            strain_displacement_matrix_w.multiply_by_matrix(&element_local_displacements)?;
        let shear_z = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_w)[0] * shear_modulus * self.area * self.shape_factor;
        forces_components.push(ForceComponent::ForceZ);
        forces_values.push(shear_z);

        let strain_displacement_matrix_thu =
            BeamAuxFunctions::strain_displacement_matrix_thu(
                self.node_1_number, self.node_2_number, r, tolerance, nodes);
        let strains_matrix_thu =
            strain_displacement_matrix_thu.multiply_by_matrix(&element_local_displacements)?;
        let moment_x = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_thu)[0] * shear_modulus * self.it;
        forces_components.push(ForceComponent::MomentX);
        forces_values.push(moment_x);

        let strain_displacement_matrix_thv =
            BeamAuxFunctions::strain_displacement_matrix_thv(
                self.node_1_number, self.node_2_number, r, tolerance, nodes);
        let strains_matrix_thv =
            strain_displacement_matrix_thv.multiply_by_matrix(&element_local_displacements)?;
        let moment_y = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_thv)[0] * self.young_modulus * self.i22;
        forces_components.push(ForceComponent::MomentY);
        forces_values.push(moment_y);

        let strain_displacement_matrix_thw =
            BeamAuxFunctions::strain_displacement_matrix_thw(
                self.node_1_number, self.node_2_number, r, tolerance, nodes);
        let strains_matrix_thw =
            strain_displacement_matrix_thw.multiply_by_matrix(&element_local_displacements)?;
        let moment_z = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_thw)[0] * self.young_modulus * self.i11;
        forces_components.push(ForceComponent::MomentZ);
        forces_values.push(moment_z);

        let element_forces = ElementForces::create(forces_values,
            forces_components);

        let element_analysis_data = ElementAnalysisData::create(
            None, None, Some(element_forces));
        Ok(element_analysis_data)

    }


    fn extract_nodes_numbers(&self) -> Vec<T>
    {
        let mut numbers = Vec::new();
        let node_1_number = self.node_1_number;
        let node_2_number = self.node_2_number;
        numbers.push(node_1_number);
        numbers.push(node_2_number);
        numbers
    }


    fn extract_fe_properties(&self) -> Vec<V>
    {
        let mut properties = Vec::new();
        properties.push(self.young_modulus);
        properties.push(self.area);
        properties
    }
}
