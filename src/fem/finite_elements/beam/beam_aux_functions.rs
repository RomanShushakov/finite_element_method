use std::ops::{Add, Sub, Div, Rem, SubAssign, Mul, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::collections::HashMap;
use std::f32::consts::PI;

use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::functions::matrix_element_value_extractor;

use extended_matrix_float::MyFloatTrait;

use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData, GlobalDOFParameter};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::beam::consts::{BEAM2N1IPT_NODES_NUMBER, BEAM_NODE_DOF};
use crate::fem::finite_elements::functions::compare_with_tolerance;


pub struct BeamAuxFunctions<T, V>(T, V);


impl<T, V> BeamAuxFunctions<T, V>
    where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
             From<u8> + Ord + 'static,
          V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + From<f32> +
             Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
             MyFloatTrait + PartialOrd + MyFloatTrait<Other = V> + 'static,
{
    pub fn length(node_1_number: T, node_2_number: T, nodes: &HashMap<T, FENode<V>>) -> V
    {
        let node_1 = nodes.get(&node_1_number).unwrap();
        let node_2 = nodes.get(&node_2_number).unwrap();

        (
            (node_1.copy_x() - node_2.copy_x()).my_powi(2) +
            (node_1.copy_y() - node_2.copy_y()).my_powi(2) +
            (node_1.copy_z() - node_2.copy_z()).my_powi(2)
        ).my_sqrt()
    }


    pub fn nodes_number() -> T
    {
        let mut n = T::from(0u8);
        (0..BEAM2N1IPT_NODES_NUMBER).for_each(|_| n += T::from(1u8));
        n
    }


    pub fn node_dof() -> T
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
            T::from(1u8), vec![a_x, a_y, a_z], tolerance)?;

        let b_x = line_b[0];
        let b_y = line_b[1];
        let b_z = line_b[2];

        let norm = V::from(1f32) / (b_x.my_powi(2) + b_y.my_powi(2) + b_z.my_powi(2));

        let mut coeff_matrix = ExtendedMatrix::create(T::from(3u8),
            T::from(3u8), vec![
                V::from(-1f32) * b_z * b_z - b_y * b_y, b_x * b_y, b_x * b_z,
                b_y * b_x, V::from(-1f32) * b_x * b_x - b_z * b_z,	b_y * b_z,
                b_z * b_x,	b_z * b_y, V::from(-1f32) * b_y * b_y - b_x * b_x,
            ], tolerance)?;

        coeff_matrix.multiply_by_number(norm);

        let components_of_line_a_perpendicular_to_line_b_matrix = coeff_matrix
            .multiply_by_matrix(&a)?;

        let a_perpendicular_to_b_x = matrix_element_value_extractor(
            T::from(0u8), T::from(0u8), &components_of_line_a_perpendicular_to_line_b_matrix)?;

        let a_perpendicular_to_b_y = matrix_element_value_extractor(
            T::from(1u8), T::from(0u8), &components_of_line_a_perpendicular_to_line_b_matrix)?;

        let a_perpendicular_to_b_z = matrix_element_value_extractor(
            T::from(2u8), T::from(0u8), &components_of_line_a_perpendicular_to_line_b_matrix)?;

        Ok([a_perpendicular_to_b_x, a_perpendicular_to_b_y, a_perpendicular_to_b_z])
    }


    pub fn rotation_matrix(node_1_number: T, node_2_number: T, local_axis_1_direction: &[V; 3],
        angle: V, tolerance: V, nodes: &HashMap<T, FENode<V>>)
        -> Result<ExtendedMatrix<T, V>, String>
    {
        let node_1 = nodes.get(&node_1_number).unwrap();
        let node_2 = nodes.get(&node_2_number).unwrap();

        let x = node_2.copy_x() - node_1.copy_x();
        let y = node_2.copy_y() - node_1.copy_y();
        let z = node_2.copy_z() - node_1.copy_z();

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
            3, 3,
            vec![q_11, q_12, q_13, q_21, q_22, q_23, q_31, q_32, q_33],
            tolerance)?;

        let projection_of_beam_section_orientation = ExtendedMatrix::create(
            3, 1,
            components_projection_of_beam_section_orientation_vector.to_vec(),
            tolerance)?;

        let transformed_projection_of_beam_section_orientation =
            interim_rotation_matrix.multiply_by_matrix(
                &projection_of_beam_section_orientation)?;

        let transformed_projection_of_beam_section_orientation_x = 
            matrix_element_value_extractor(0, 0, &transformed_projection_of_beam_section_orientation)?;

        let transformed_projection_of_beam_section_orientation_y =
            matrix_element_value_extractor(1, 0, &transformed_projection_of_beam_section_orientation)?;

        let transformed_projection_of_beam_section_orientation_z =
            matrix_element_value_extractor(2, 0, &transformed_projection_of_beam_section_orientation)?;

        let angle_between_beam_section_local_axis_1_direction_and_axis_t =
            (transformed_projection_of_beam_section_orientation_z /
            (transformed_projection_of_beam_section_orientation_x.my_powi(2) +
            transformed_projection_of_beam_section_orientation_y.my_powi(2) +
            transformed_projection_of_beam_section_orientation_z.my_powi(2))
                .my_sqrt()
            ).my_acos();

        let total_angle = angle +
            angle_between_beam_section_local_axis_1_direction_and_axis_t;

        let c_x = compare_with_tolerance(x / length, tolerance);
        let c_y = compare_with_tolerance(y / length, tolerance);
        let c_z = compare_with_tolerance(z / length, tolerance);
        let c_xz = compare_with_tolerance(
            (c_x.my_powi(2) + c_z.my_powi(2)).my_sqrt(), tolerance);

        let c =  compare_with_tolerance(total_angle.my_cos(), tolerance);
        let s = compare_with_tolerance(total_angle.my_sin(), tolerance);

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


    fn inverse_jacobian(node_1_number: T, node_2_number: T, r: V, nodes: &HashMap<T, FENode<V>>)
        -> V
    {
        V::from(1f32) / BeamAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
    }


    pub fn determinant_of_jacobian(node_1_number: T, node_2_number: T, r: V,
        nodes: &HashMap<T, FENode<V>>) -> V
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


    pub fn strain_displacement_matrix_u(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let elements = vec![
                BeamAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
                BeamAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
            ];
        let mut matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            elements, tolerance)?;
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        matrix.multiply_by_number(inverse_jacobian);
        Ok(matrix)
    }


    pub fn strain_displacement_matrix_v(node_1_number: T, node_2_number: T, r: V, tolerance: V,
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
            lhs_elements, tolerance)?;
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
            rhs_elements, tolerance)?;

        let matrix = lhs_matrix.subtract_matrix(&rhs_matrix)?;

        Ok(matrix)
    }


    pub fn strain_displacement_matrix_w(node_1_number: T, node_2_number: T, r: V, tolerance: V,
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
            lhs_elements, tolerance)?;
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
            rhs_elements, tolerance)?;

        let matrix = lhs_matrix.subtract_matrix(&rhs_matrix)?;

        Ok(matrix)
    }


    pub fn strain_displacement_matrix_thu(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let elements = vec![
                V::from(0f32), V::from(0f32), V::from(0f32),
                BeamAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
                BeamAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32), V::from(0f32),
            ];
        let mut matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            elements, tolerance)?;
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        matrix.multiply_by_number(inverse_jacobian);
        Ok(matrix)
    }


    pub fn strain_displacement_matrix_thv(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let elements = vec![
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), BeamAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), BeamAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32),
            ];
        let mut matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            elements, tolerance)?;
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        matrix.multiply_by_number(inverse_jacobian);
        Ok(matrix)
    }


    pub fn strain_displacement_matrix_thw(node_1_number: T, node_2_number: T, r: V, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let elements = vec![
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), BeamAuxFunctions::<T, V>::dh1_dr(r),
                V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), BeamAuxFunctions::<T, V>::dh2_dr(r),
            ];
        let mut matrix = ExtendedMatrix::create(T::from(1u8),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            elements, tolerance)?;
        let inverse_jacobian = BeamAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
            r, nodes);
        matrix.multiply_by_number(inverse_jacobian);
        Ok(matrix)
    }


    pub fn local_stiffness_matrix(node_1_number: T, node_2_number: T, young_modulus: V,
        poisson_ratio: V, area: V, shear_factor: V, it: V, i11: V, i22: V,
        alpha: V, r: V, tolerance: V, nodes: &HashMap<T, FENode<V>>)
        -> Result<ExtendedMatrix<T, V>, String>
    {
        let mut lhs_matrix_u = BeamAuxFunctions::strain_displacement_matrix_u(
            node_1_number, node_2_number, r, tolerance, nodes)?;
        lhs_matrix_u.transpose();
        let rhs_matrix_u = BeamAuxFunctions::strain_displacement_matrix_u(
            node_1_number, node_2_number, r, tolerance, nodes)?;
        let mut matrix_u = lhs_matrix_u.multiply_by_matrix(&rhs_matrix_u)
            .map_err(|e|
                format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let coeff_u = young_modulus * area;
        matrix_u.multiply_by_number(coeff_u);
        matrix_u.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

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
        let coeff_v = shear_modulus * area * shear_factor;
        matrix_v.multiply_by_number(coeff_v);
        matrix_v.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

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
        let coeff_w = shear_modulus * area * shear_factor;
        matrix_w.multiply_by_number(coeff_w);
        matrix_w.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

        let mut lhs_matrix_thu = BeamAuxFunctions::strain_displacement_matrix_thu(
            node_1_number, node_2_number, r, tolerance, nodes)?;
        lhs_matrix_thu.transpose();
        let rhs_matrix_thu = BeamAuxFunctions::strain_displacement_matrix_thu(
            node_1_number, node_2_number, r, tolerance, nodes)?;
        let mut matrix_thu = lhs_matrix_thu.multiply_by_matrix(
            &rhs_matrix_thu)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let coeff_thu = shear_modulus * it;
        matrix_thu.multiply_by_number(coeff_thu);
        matrix_thu.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

        let mut lhs_matrix_thv = BeamAuxFunctions::strain_displacement_matrix_thv(
            node_1_number, node_2_number, r, tolerance, nodes)?;
        lhs_matrix_thv.transpose();
        let rhs_matrix_thv = BeamAuxFunctions::strain_displacement_matrix_thv(
            node_1_number, node_2_number, r, tolerance, nodes)?;
        let mut matrix_thv = lhs_matrix_thv.multiply_by_matrix(
            &rhs_matrix_thv)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let coeff_thv = young_modulus * i22;
        matrix_thv.multiply_by_number(coeff_thv);
        matrix_thv.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

        let mut lhs_matrix_thw = BeamAuxFunctions::strain_displacement_matrix_thw(
            node_1_number, node_2_number, r, tolerance, nodes)?;
        lhs_matrix_thw.transpose();
        let rhs_matrix_thw = BeamAuxFunctions::strain_displacement_matrix_thw(
            node_1_number, node_2_number, r, tolerance, nodes)?;
        let mut matrix_thw = lhs_matrix_thw.multiply_by_matrix(
            &rhs_matrix_thw)
            .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be calculated! \
                Reason: {}", e))?;
        let coeff_thw = young_modulus * i11;
        matrix_thw.multiply_by_number(coeff_thw);
        matrix_thw.multiply_by_number(BeamAuxFunctions::determinant_of_jacobian(
            node_1_number, node_2_number, r, nodes) * alpha);

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


    pub fn compose_node_dof_parameters<'a>(node_number: T)
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
                let value = matrix_element_value_extractor(row, column, &column_matrix)?;
                values.push(value);
                column += T::from(1u8);
            }
            row += T::from(1u8);
        }
        Ok(values)
    }


    pub fn find_principal_moments_of_inertia(i11_init: V, i22_init: V, i12_init: V) -> (V, V, V)
    {
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
