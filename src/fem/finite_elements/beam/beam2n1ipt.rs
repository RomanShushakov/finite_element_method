// use std::hash::Hash;
// use std::fmt::Debug;
// use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};
//
// use extended_matrix::basic_matrix::basic_matrix::MatrixElementPosition;
// use extended_matrix::extended_matrix::ExtendedMatrix;
// use extended_matrix::functions::extract_element_value;
//
// use crate::fem::finite_elements::finite_element::FiniteElementTrait;
// use crate::fem::finite_elements::fe_node::FENode;
//
// use crate::fem::global_analysis::fe_stiffness::{StiffnessGroup, StiffnessType};
// use crate::fem::global_analysis::fe_dof_parameter_data::{GlobalDOFParameter, DOFParameterData};
// use crate::fem::global_analysis::fe_global_analysis_result::Displacements;
//
// use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
// use crate::fem::element_analysis::fe_element_analysis_result::ElementAnalysisData;
// use crate::fem::element_analysis::fe_stress_strain_components::StressStrainComponent;
//
// use crate::fem::finite_elements::functions::compare_with_tolerance;
//
// use crate::my_float::MyFloatTrait;
//
//
// use crate::fem::finite_elements::truss::consts::
// {
//     TRUSS_NODE_DOF, TRUSS2N1IP_NODES_NUMBER, TRUSS_STRESS_STRAIN_COMPONENTS_NUMBERS,
//     POINTS_NUMBER_FOR_TAPERED_TRUSS,
// };
// use std::collections::HashMap;
//
//
// struct TrussAuxFunctions<T, V>(T, V);
//
//
// impl<T, V> TrussAuxFunctions<T, V>
//     where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
//              Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
//              From<u8> + 'static,
//           V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + From<f32> +
//              Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
//              MyFloatTrait + PartialOrd + 'static,
// {
//     fn length(node_1_number: T, node_2_number: T, nodes: &HashMap<T, FENode<V>>) -> V
//     {
//         let node_1 = nodes.get(&node_1_number).unwrap();
//         let node_2 = nodes.get(&node_2_number).unwrap();
//
//         ((node_1.x() - node_2.x()).my_powi(2) + (node_1.y() - node_2.y()).my_powi(2) +
//         (node_1.z() - node_2.z()).my_powi(2)).my_sqrt()
//     }
//
//
//     fn nodes_number() -> T
//     {
//         let mut n = T::from(0u8);
//         (0..TRUSS2N1IP_NODES_NUMBER).for_each(|_| n += T::from(1u8));
//         n
//     }
//
//
//     fn node_dof() -> T
//     {
//         let mut n = T::from(0u8);
//         (0..TRUSS_NODE_DOF).for_each(|_| n += T::from(1u8));
//         n
//     }
//
//
//     fn rotation_matrix(node_1_number: T, node_2_number: T, tolerance: V,
//         nodes: &HashMap<T, FENode<V>>) -> ExtendedMatrix<T, V>
//     {
//         let node_1 = nodes.get(&node_1_number).unwrap();
//         let node_2 = nodes.get(&node_2_number).unwrap();
//
//         let x = node_2.x() - node_1.x();
//         let y = node_2.y() - node_1.y();
//         let z = node_2.z() - node_1.z();
//
//         let length = TrussAuxFunctions::<T, V>::length(node_1_number, node_2_number, nodes);
//
//         let (u, v, w) = (length, V::from(0f32), V::from(0f32));
//         let alpha = ((x * u + y * v + z * w) / (length * length)).my_acos();
//         let (rotation_axis_coord_x, mut rotation_axis_coord_y,
//             mut rotation_axis_coord_z) = (V::from(0f32), V::from(0f32), V::from(0f32));
//         if x != V::from(0f32) && y == V::from(0f32) && z == V::from(0f32)
//         {
//             rotation_axis_coord_z = x;
//         }
//         else
//         {
//             rotation_axis_coord_y = z * length;
//             rotation_axis_coord_z = y * V::from(-1f32) * length;
//         }
//         let norm = V::from(1f32) / (rotation_axis_coord_x.my_powi(2) +
//             rotation_axis_coord_y.my_powi(2) + rotation_axis_coord_z.my_powi(2)).my_sqrt();
//         let (x_n, y_n, z_n) = (rotation_axis_coord_x * norm,
//             rotation_axis_coord_y * norm, rotation_axis_coord_z * norm);
//         let (c, s) = (alpha.my_cos(), alpha.my_sin());
//         let t = V::from(1f32) - c;
//         let q_11 = compare_with_tolerance(t * x_n * x_n + c, tolerance);
//         let q_12 = compare_with_tolerance(t * x_n * y_n - z_n * s, tolerance);
//         let q_13 = compare_with_tolerance(t * x_n * z_n + y_n * s, tolerance);
//         let q_21 = compare_with_tolerance(t * x_n * y_n + z_n * s, tolerance);
//         let q_22 = compare_with_tolerance(t * y_n * y_n + c, tolerance);
//         let q_23 = compare_with_tolerance(t * y_n * z_n - x_n * s, tolerance);
//         let q_31 = compare_with_tolerance(t * x_n * z_n - y_n * s, tolerance);
//         let q_32 = compare_with_tolerance(t * y_n * z_n + x_n * s, tolerance);
//         let q_33 = compare_with_tolerance(t * z_n * z_n + c, tolerance);
//
//         ExtendedMatrix::create(
//             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//             vec![
//                 [q_11, q_12, q_13], [V::from(0f32); TRUSS_NODE_DOF],
//                 [q_21, q_22, q_23], [V::from(0f32); TRUSS_NODE_DOF],
//                 [q_31, q_32, q_33], [V::from(0f32); TRUSS_NODE_DOF],
//                 [V::from(0f32); TRUSS_NODE_DOF], [q_11, q_12, q_13],
//                 [V::from(0f32); TRUSS_NODE_DOF], [q_21, q_22, q_23],
//                 [V::from(0f32); TRUSS_NODE_DOF], [q_31, q_32, q_33],
//             ].concat(),
//             tolerance,
//         )
//     }
//
//
//     fn power_func_x(a: V, x: V, n: i32) -> V
//     {
//         (0..n).fold(a, |acc, _| acc * x)
//     }
//
//
//     fn derivative_x(f: fn(V, V, i32) -> V, a: V, x: V, n: i32) -> V
//     {
//         let mut converted_n = V::from(0f32);
//         (0..n).for_each(|_| converted_n += V::from(1f32));
//
//         f(a * converted_n, x, n - 1)
//     }
//
//
//     fn dx_dr(x_1: V, x_2: V, r: V) -> V
//     {
//         TrussAuxFunctions::<T, V>::derivative_x(
//             TrussAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.5f32), V::from(0f32),
//             0) -
//         TrussAuxFunctions::<T, V>::derivative_x(
//             TrussAuxFunctions::<T, V>::power_func_x, x_1 * V::from(0.5f32), r, 1) +
//         TrussAuxFunctions::<T, V>::derivative_x(
//             TrussAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.5f32), V::from(0f32),
//             0) +
//         TrussAuxFunctions::<T, V>::derivative_x(
//             TrussAuxFunctions::<T, V>::power_func_x, x_2 * V::from(0.5f32), r, 1)
//     }
//
//
//     fn jacobian(node_1_number: T, node_2_number: T, r: V, nodes: &HashMap<T, FENode<V>>) -> V
//     {
//         let length = TrussAuxFunctions::length(node_1_number, node_2_number, nodes);
//         let x_1 = V::from(-1f32) * length / V::from(2f32);
//         let x_2 = length / V::from(2f32);
//         TrussAuxFunctions::<T, V>::dx_dr(x_1, x_2, r)
//     }
//
//
//     fn inverse_jacobian(node_1_number: T, node_2_number: T, r: V, nodes: &HashMap<T, FENode<V>>) -> V
//     {
//         V::from(1f32) / TrussAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
//     }
//
//
//     fn determinant_of_jacobian(node_1_number: T, node_2_number: T, r: V,
//         nodes: &HashMap<T, FENode<V>>) -> V
//     {
//         TrussAuxFunctions::jacobian(node_1_number, node_2_number, r, nodes)
//     }
//
//
//     fn dh1_dr(r: V) -> V
//     {
//         TrussAuxFunctions::<T, V>::derivative_x(
//             TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) -
//         TrussAuxFunctions::<T, V>::derivative_x(
//             TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), r, 1)
//     }
//
//
//     fn dh2_dr(r: V) -> V
//     {
//         TrussAuxFunctions::<T, V>::derivative_x(
//             TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), V::from(0f32), 0) +
//         TrussAuxFunctions::<T, V>::derivative_x(
//             TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5f32), r, 1)
//     }
//
//
//     fn strain_displacement_matrix(node_1_number: T, node_2_number: T, r: V, tolerance: V,
//         nodes: &HashMap<T, FENode<V>>) -> ExtendedMatrix<T, V>
//     {
//         let elements = vec![TrussAuxFunctions::<T, V>::dh1_dr(r), V::from(0f32),
//             V::from(0f32), TrussAuxFunctions::<T, V>::dh2_dr(r), V::from(0f32), V::from(0f32)];
//         let mut matrix = ExtendedMatrix::create(T::from(1u8),
//             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//             elements, tolerance);
//         let inverse_jacobian = TrussAuxFunctions::inverse_jacobian(node_1_number, node_2_number,
//             r, nodes);
//         matrix.multiply_by_number(inverse_jacobian);
//         matrix
//     }
//
//
//     fn area(area_1: V, area_2: Option<V>, r: V) -> V
//     {
//         if let Some(area_2) = area_2
//         {
//             // V::from((area_1.into().sqrt() + (area_2.into().sqrt() - area_1.into().sqrt()) *
//             //     (r.into() + 1.0) / 2.0).powi(2))
//             (area_2 - area_1) / V::from(2f32) * r + area_1 -
//                 (area_2 - area_1) / V::from(2f32) * V::from(-1f32)
//         }
//         else
//         {
//             area_1
//         }
//     }
//
//
//     fn local_stiffness_matrix(node_1_number: T, node_2_number: T, young_modulus: V, area_1: V,
//         area_2: Option<V>, alpha: V, r: V, local_stiffness_matrix: &ExtendedMatrix<T, V>,
//         tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<ExtendedMatrix<T, V>, String>
//     {
//         let current_area = TrussAuxFunctions::<T, V>::area(area_1, area_2, r);
//
//         let mut lhs_matrix = TrussAuxFunctions::strain_displacement_matrix(
//             node_1_number, node_2_number, r, tolerance, nodes);
//
//         lhs_matrix.transpose();
//
//         lhs_matrix.multiply_by_number(young_modulus * current_area);
//
//         let rhs_matrix = TrussAuxFunctions::strain_displacement_matrix(
//                 node_1_number, node_2_number, r, tolerance, nodes);
//
//         return match lhs_matrix.multiply_by_matrix(&rhs_matrix)
//         {
//             Ok(mut matrix) =>
//                 {
//                     matrix.multiply_by_number(TrussAuxFunctions::determinant_of_jacobian(
//                         node_1_number, node_2_number, r, nodes) * alpha);
//
//                     match local_stiffness_matrix.add_matrix(&matrix)
//                     {
//                         Ok(matrix) => Ok(matrix),
//                         Err(e) =>
//                             Err(format!("Truss2n2ip: Local stiffness matrix cannot be \
//                                 calculated! Reason: {}", e)),
//                     }
//                 },
//             Err(e) => Err(format!("Truss2n2ip: Local stiffness matrix cannot be \
//                                 calculated! Reason: {}", e)),
//         }
//     }
//
//
//     fn compose_node_dof_parameters<'a>(node_number: T)
//         -> Result<Vec<DOFParameterData<T>>, &'a str>
//     {
//         let mut node_dof_parameters = Vec::new();
//         for dof in 0..TRUSS_NODE_DOF
//         {
//             let dof_parameter =
//                 GlobalDOFParameter::iterator().nth(dof)
//                     .ok_or("Truss2n2ip: Could not find dof parameter!")?;
//             let dof_parameter = DOFParameterData { node_number,
//                 dof_parameter: *dof_parameter
//             };
//             node_dof_parameters.push(dof_parameter);
//         }
//         Ok(node_dof_parameters)
//     }
//
//
//     fn extract_column_matrix_values(column_matrix: &ExtendedMatrix<T, V>) -> Vec<V>
//     {
//         let mut values = Vec::new();
//         let shape = column_matrix.get_shape();
//         let all_values = column_matrix.extract_all_elements_values();
//
//         let mut row = T::from(0u8);
//         while row < shape.0
//         {
//             let mut column = T::from(0u8);
//             while column < shape.1
//             {
//                 let value = extract_element_value(row, column, &all_values);
//                 values.push(value);
//                 column += T::from(1u8);
//             }
//             row += T::from(1u8);
//         }
//         values
//     }
// }
//
//
// struct IntegrationPoint<V>
// {
//     r: V,
//     weight: V,
// }
//
//
// struct State<T, V>
// {
//     rotation_matrix: ExtendedMatrix<T, V>,
//     integration_points: Vec<IntegrationPoint<V>>,
//     local_stiffness_matrix: ExtendedMatrix<T, V>,
//     nodes_dof_parameters_global: Vec<DOFParameterData<T>>,
// }
//
//
// impl<T, V> State<T, V>
// {
//     fn create(rotation_matrix: ExtendedMatrix<T, V>, integration_points: Vec<IntegrationPoint<V>>,
//         local_stiffness_matrix: ExtendedMatrix<T, V>,
//         nodes_dof_parameters_global: Vec<DOFParameterData<T>>) -> Self
//     {
//         State { rotation_matrix, integration_points, local_stiffness_matrix,
//         nodes_dof_parameters_global }
//     }
// }
//
//
// pub struct Truss2n1ip<T, V>
// {
//     node_1_number: T,
//     node_2_number: T,
//     young_modulus: V,
//     area: V,
//     area_2: Option<V>,
//     state: State<T, V>,
// }
//
//
// impl<T, V> Truss2n1ip<T, V>
//     where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
//              Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
//              From<u8> + 'static,
//           V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + From<f32> + Add<Output = V> +
//              Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
//              MyFloatTrait + PartialOrd + 'static
// {
//     pub fn create(node_1_number: T, node_2_number: T, young_modulus: V, area: V, area_2: Option<V>,
//         tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<Self, String>
//     {
//         let integration_point_1 = IntegrationPoint {
//             r: V::from(0f32), weight: V::from(2f32) };
//
//         let rotation_matrix = TrussAuxFunctions::<T, V>::rotation_matrix(
//             node_1_number, node_2_number, tolerance, nodes);
//
//         let integration_points = vec![integration_point_1];
//
//         let mut local_stiffness_matrix = ExtendedMatrix::create(
//             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//             vec![V::from(0f32); (TRUSS2N1IP_NODES_NUMBER * TRUSS_NODE_DOF).pow(2)], tolerance);
//
//         for integration_point in &integration_points
//         {
//             let matrix = TrussAuxFunctions::local_stiffness_matrix(
//                 node_1_number, node_2_number, young_modulus, area, area_2,
//                 integration_point.weight, integration_point.r, &local_stiffness_matrix,
//                 tolerance, nodes)?;
//             local_stiffness_matrix = matrix;
//         }
//
//         let mut nodes_dof_parameters =
//             TrussAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;
//
//         let node_2_dof_parameters =
//             TrussAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;
//
//         nodes_dof_parameters.extend(node_2_dof_parameters);
//
//         let state = State::create(rotation_matrix, integration_points,
//             local_stiffness_matrix, nodes_dof_parameters);
//
//         Ok(Truss2n1ip { node_1_number, node_2_number, young_modulus, area, area_2, state })
//     }
//
//
//     fn extract_local_displacements(&self, global_displacements: &Displacements<T, V>, tolerance: V)
//         -> Result<ExtendedMatrix<T, V>, String>
//     {
//         let mut element_global_displacements_values = Vec::new();
//         for lhs_dof_parameter_data in &self.state.nodes_dof_parameters_global
//         {
//             if let Some(position) = global_displacements.dof_parameters_data
//                 .iter()
//                 .position(|rhs_dof_parameter_data|
//                     rhs_dof_parameter_data == lhs_dof_parameter_data)
//             {
//                 let displacement_value = global_displacements.displacements_values[position];
//                 element_global_displacements_values.push(displacement_value);
//             }
//             else
//             {
//                 element_global_displacements_values.push(V::from(0f32));
//             }
//         }
//
//         let mut m = 0usize;
//         let mut rows_number = T::from(0u8);
//         while m < self.state.nodes_dof_parameters_global.len()
//         {
//             m += 1usize;
//             rows_number += T::from(1u8);
//         }
//
//         let element_global_displacements = ExtendedMatrix::create(rows_number,
//             T::from(1u8), element_global_displacements_values,
//             tolerance);
//
//         let element_local_displacements =
//             self.state.rotation_matrix.multiply_by_matrix(&element_global_displacements)?;
//         Ok(element_local_displacements)
//     }
// }
//
//
// impl<T, V> FiniteElementTrait<T, V> for Truss2n1ip<T, V>
//     where T: Copy + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> +
//              Mul<Output = T> + Eq + Hash + Debug + SubAssign + PartialOrd + AddAssign +
//              From<u8> + 'static,
//           V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
//              Into<f64> + SubAssign + AddAssign + MulAssign + PartialEq + Debug +
//              MyFloatTrait + PartialOrd + From<f32> + 'static,
// {
//     fn update(&mut self, nodes_numbers: Vec<T>, properties: Vec<V>, tolerance: V,
//         nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
//     {
//         let node_1_number = nodes_numbers[0];
//
//         let node_2_number = nodes_numbers[1];
//
//         let young_modulus = properties[0];
//
//         let area = properties[1];
//
//         let area_2 =
//             if properties.len() == 3 { Some(properties[2]) } else { None };
//
//         let rotation_matrix = TrussAuxFunctions::rotation_matrix(
//             node_1_number, node_2_number, tolerance, nodes);
//
//         let mut local_stiffness_matrix = ExtendedMatrix::create(
//             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//             vec![V::from(0f32); (TRUSS2N1IP_NODES_NUMBER * TRUSS_NODE_DOF).pow(2)], tolerance);
//
//         for integration_point in &self.state.integration_points
//         {
//             let matrix = TrussAuxFunctions::local_stiffness_matrix(
//                 node_1_number, node_2_number, young_modulus, area, area_2,
//                 integration_point.weight, integration_point.r, &local_stiffness_matrix,
//                 tolerance, nodes)?;
//             local_stiffness_matrix = matrix;
//         }
//
//         let mut nodes_dof_parameters =
//             TrussAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;
//
//         let node_2_dof_parameters =
//             TrussAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;
//
//         nodes_dof_parameters.extend(node_2_dof_parameters);
//
//         self.node_1_number = node_1_number;
//         self.node_2_number = node_2_number;
//         self.young_modulus = young_modulus;
//         self.area = area;
//         self.area_2 = area_2;
//         self.state.rotation_matrix = rotation_matrix;
//         self.state.local_stiffness_matrix = local_stiffness_matrix;
//         self.state.nodes_dof_parameters_global = nodes_dof_parameters;
//
//         Ok(())
//     }
//
//
//     fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
//     {
//         let mut interim_matrix = self.state.rotation_matrix.clone();
//         interim_matrix.transpose();
//         if let Ok(matrix) =
//         interim_matrix.multiply_by_matrix(&self.state.local_stiffness_matrix)
//         {
//             if let Ok(matrix) =
//             matrix.multiply_by_matrix(&self.state.rotation_matrix)
//             {
//                 return Ok(matrix);
//             }
//         }
//         Err("Truss2n2ip: Stiffness matrix cannot be extracted!")
//     }
//
//
//     fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>
//     {
//         let (rows_number, columns_number) =
//             (TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//              TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof());
//
//         let mut positions_1_1 = Vec::new();
//         let mut positions_1_2 = Vec::new();
//         let mut positions_2_1 = Vec::new();
//         let mut positions_2_2 = Vec::new();
//
//         let mut i = T::from(0u8);
//         while i < rows_number * columns_number
//         {
//             let position = MatrixElementPosition::create(
//                 i / columns_number, i % columns_number);
//
//             let row = i / columns_number;
//             let column = i % columns_number;
//
//             if row < TrussAuxFunctions::<T, V>::node_dof() && column < TrussAuxFunctions::<T, V>::node_dof()
//             {
//                 positions_1_1.push(position);
//             }
//             else if row < TrussAuxFunctions::<T, V>::node_dof() && column >= TrussAuxFunctions::<T, V>::node_dof()
//             {
//                 positions_1_2.push(position);
//             }
//             else if row >= TrussAuxFunctions::<T, V>::node_dof() && column < TrussAuxFunctions::<T, V>::node_dof()
//             {
//                 positions_2_1.push(position);
//             }
//             else
//             {
//                 positions_2_2.push(position);
//             }
//             i += T::from(1u8);
//         }
//
//         vec![StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
//                 number_2: self.node_1_number, positions: positions_1_1, },
//              StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
//                 number_2: self.node_2_number, positions: positions_1_2, },
//              StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
//                 number_2: self.node_1_number, positions: positions_2_1 },
//              StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
//                 number_2: self.node_2_number, positions: positions_2_2 },
//         ]
//     }
//
//
//     fn node_belong_element(&self, node_number: T) -> bool
//     {
//         self.node_1_number == node_number || self.node_2_number == node_number
//     }
//
//
//     fn refresh(&mut self, tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
//     {
//         let rotation_matrix = TrussAuxFunctions::rotation_matrix(
//             self.node_1_number, self.node_2_number, tolerance, nodes);
//
//         let mut local_stiffness_matrix = ExtendedMatrix::create(
//             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
//             vec![V::from(0f32); (TRUSS2N1IP_NODES_NUMBER * TRUSS_NODE_DOF).pow(2)], tolerance);
//
//         for integration_point in self.state.integration_points.iter()
//         {
//             let matrix = TrussAuxFunctions::local_stiffness_matrix(
//                 self.node_1_number, self.node_2_number, self.young_modulus, self.area,
//                 self.area_2, integration_point.weight, integration_point.r,
//                 &local_stiffness_matrix, tolerance, nodes)?;
//             local_stiffness_matrix = matrix;
//         }
//
//         self.state.rotation_matrix = rotation_matrix;
//         self.state.local_stiffness_matrix = local_stiffness_matrix;
//         Ok(())
//     }
//
//
//     fn nodes_numbers_same(&self, nodes_numbers: Vec<T>) -> bool
//     {
//         (nodes_numbers[0] == self.node_1_number && nodes_numbers[1] == self.node_2_number) ||
//         (nodes_numbers[0] == self.node_2_number && nodes_numbers[1] == self.node_1_number)
//     }
//
//
//     fn extract_element_analysis_data(&self, global_displacements: &Displacements<T, V>,
//         tolerance: V, nodes: &HashMap<T, FENode<V>>, number: T)
//         -> Result<ElementAnalysisData<T, V>, String>
//     {
//         let element_local_displacements =
//             self.extract_local_displacements(global_displacements, tolerance)?;
//         if self.area_2.is_some()
//         {
//             let mut forces_components = Vec::new();
//             let mut forces_values = Vec::new();
//             let mut axial_force = V::from(0f32);
//             for ip in self.state.integration_points.iter()
//             {
//                 let strain_displacement_matrix =
//                     TrussAuxFunctions::strain_displacement_matrix(self.node_1_number,
//                         self.node_2_number, ip.r, tolerance, nodes);
//
//                 let strains_matrix =
//                     strain_displacement_matrix.multiply_by_matrix(
//                         &element_local_displacements)?;
//
//                 let stresses_matrix =
//                     {
//                         let mut matrix = strains_matrix.clone();
//                         matrix.multiply_by_number(self.young_modulus);
//                         matrix
//                     };
//                 let stresses_values =
//                     TrussAuxFunctions::extract_column_matrix_values(&stresses_matrix);
//                 for stress in &stresses_values
//                 {
//                     let area = TrussAuxFunctions::<T, V>::area(
//                         self.area, self.area_2, ip.r);
//                     axial_force += *stress * area;
//                 }
//             }
//
//             let mut coeff = V::from(0f32);
//             (0..self.state.integration_points.len()).for_each(|_| coeff += V::from(1f32));
//             let integral_axial_force = axial_force / coeff;
//
//             forces_components.push(ForceComponent::Axial);
//             forces_values.push(integral_axial_force);
//             let mut strains_components = Vec::new();
//             let mut strains_values = Vec::new();
//             let mut stresses_components = Vec::new();
//             let mut stresses_values = Vec::new();
//             for point_number in 0..POINTS_NUMBER_FOR_TAPERED_TRUSS
//             {
//                 for component_number in &TRUSS_STRESS_STRAIN_COMPONENTS_NUMBERS
//                 {
//                     let stress_strain_component =
//                         StressStrainComponent::iterator().nth(*component_number as usize)
//                             .ok_or(format!("Truss2n2ip: Stress strain component number {} \
//                                 could not be extracted", component_number))?;
//                     strains_components.push(*stress_strain_component);
//                     stresses_components.push(*stress_strain_component);
//                 }
//
//
//                 let mut coeff_1 = V::from(0f32);
//                 (0..POINTS_NUMBER_FOR_TAPERED_TRUSS - 1).for_each(|_| coeff_1 += V::from(1f32));
//                 let mut coeff_2 = V::from(0f32);
//                 (0..point_number).for_each(|_| coeff_2 += V::from(1f32));
//
//                 let r = V::from(-1f32) + V::from(2f32) / coeff_1 * coeff_2;
//
//                 let area = TrussAuxFunctions::<T, V>::area(self.area, self.area_2, r);
//                 let stress_value = integral_axial_force / area;
//                 stresses_values.push(stress_value);
//                 let strain_value = stress_value / self.young_modulus;
//                 strains_values.push(strain_value);
//             }
//             let element_analysis_data = ElementAnalysisData::create(
//                 number, strains_values, strains_components,
//                 stresses_values, stresses_components, forces_values, forces_components);
//             Ok(element_analysis_data)
//         }
//         else
//         {
//             let r = self.state.integration_points[0].r;
//             let mut strains_components = Vec::new();
//             let mut stresses_components = Vec::new();
//             let mut forces_components = Vec::new();
//             let mut forces_values = Vec::new();
//             let strain_displacement_matrix =
//                 TrussAuxFunctions::strain_displacement_matrix(
//                     self.node_1_number, self.node_2_number, r, tolerance, nodes);
//
//             let strains_matrix =
//                 strain_displacement_matrix.multiply_by_matrix(&element_local_displacements)?;
//
//             let stresses_matrix =
//                 {
//                     let mut matrix = strains_matrix.clone();
//                     matrix.multiply_by_number(self.young_modulus);
//                     matrix
//                 };
//
//             for component_number in &TRUSS_STRESS_STRAIN_COMPONENTS_NUMBERS
//             {
//                 let stress_strain_component =
//                     StressStrainComponent::iterator().nth(*component_number)
//                         .ok_or(format!("Truss2n2ip: Stress strain component number {} could \
//                             not be extracted", component_number))?;
//                 strains_components.push(*stress_strain_component);
//                 stresses_components.push(*stress_strain_component);
//             }
//             let strains_values =
//                 TrussAuxFunctions::extract_column_matrix_values(&strains_matrix);
//             let stresses_values =
//                 TrussAuxFunctions::extract_column_matrix_values(&stresses_matrix);
//             for stress in &stresses_values
//             {
//                 let area = TrussAuxFunctions::<T, V>::area(self.area, self.area_2, r);
//                 let axial_force = *stress * area;
//                 forces_components.push(ForceComponent::Axial);
//                 forces_values.push(axial_force);
//             }
//             let element_analysis_data = ElementAnalysisData::create(
//                 number, strains_values, strains_components,
//                 stresses_values, stresses_components, forces_values, forces_components);
//             Ok(element_analysis_data)
//         }
//     }
//
//
//     fn extract_nodes_numbers(&self) -> Vec<T>
//     {
//         let mut numbers = Vec::new();
//         let node_1_number = self.node_1_number;
//         let node_2_number = self.node_2_number;
//         numbers.push(node_1_number);
//         numbers.push(node_2_number);
//         numbers
//     }
//
//
//     fn extract_fe_properties(&self) -> Vec<V>
//     {
//         let mut properties = Vec::new();
//         properties.push(self.young_modulus);
//         properties.push(self.area);
//         if let Some(area) = self.area_2
//         {
//             properties.push(area);
//         }
//         properties
//     }
// }
