use std::fmt::Debug;
use std::collections::HashMap;
use std::f32::consts::PI;

use extended_matrix::
{
    FloatTrait, Vector3, VectorTrait, BasicOperationsTrait, Position, Matrix, TryIntoSquareMatrixTrait, SquareMatrix,
    Vector
};

use crate::fem::structs::{Node, NODE_DOF};
use crate::fem::methods_for_element_analysis::ElementForceComponent;
use crate::fem::quadrilateral_4n_element_functions::
{
    is_points_of_quadrilateral_on_the_same_line, is_points_of_quadrilateral_on_the_same_plane, 
    find_rotation_matrix_elements_of_quadrilateral, convex_hull_on_four_points_on_plane,
};


const PLATE_NODES_NUMBER: usize = 4;
pub const PLATE_NODE_DOF: usize = 6;


enum PlateElementDataError<V>
{
    YoungModulus(V),
    PoissonRatio(V),
    Thickness(V),
    ShearFactor(V),
    NodesOnLine(u32),
    NodesNotOnPlane(u32),
    NotConvex(u32),
}


impl<V> PlateElementDataError<V>
    where V: Debug
{
    fn compose_error_message(&self) -> String
    {
        match self
        {
            Self::YoungModulus(value) => format!("Young's modulus {value:?} is less or equal to zero!"),
            Self::PoissonRatio(value) => format!("Poisson's ratio {value:?} is less or equal to zero!"),
            Self::Thickness(value) => format!("Thickness {value:?} is less or equal to zero!"),
            Self::ShearFactor(value) => format!("Shear factor {value:?} is less or equal to zero!"),
            Self::NodesOnLine(n) => format!("Some nodes of {n} element lie on the line!"),
            Self::NodesNotOnPlane(n) => format!("Not all nodes of element {n} lie on the plane!"),
            Self::NotConvex(n) => format!("Element {n} non-convex!"),
        }
    }
}


fn check_plate_properties<V>(
    number: u32,
    young_modulus: V,
    poisson_ratio: V,
    thickness: V,
    shear_factor: V,
    node_1_number: u32,
    node_2_number: u32,
    node_3_number: u32,
    node_4_number: u32,
    nodes: &HashMap<u32, Node<V>>,
    rel_tol: V,
    abs_tol: V,
) 
    -> Result<(), String>
    where V: FloatTrait<Output = V>
{
    if young_modulus <= V::from(0f32)
    {
        return Err(PlateElementDataError::<V>::YoungModulus(young_modulus).compose_error_message());
    }
    if poisson_ratio <= V::from(0f32)
    {
        return Err(PlateElementDataError::<V>::PoissonRatio(young_modulus).compose_error_message());
    }
    if thickness <= V::from(0f32)
    {
        return Err(PlateElementDataError::<V>::Thickness(young_modulus).compose_error_message());
    }
    if shear_factor <= V::from(0f32)
    {
        return Err(PlateElementDataError::<V>::ShearFactor(shear_factor).compose_error_message());
    }
    if is_points_of_quadrilateral_on_the_same_line(
        &nodes.get(&node_1_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
        &nodes.get(&node_2_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
        &nodes.get(&node_3_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
        &nodes.get(&node_4_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
    )
    {
        return Err(PlateElementDataError::<V>::NodesOnLine(number).compose_error_message());
    }
    if !is_points_of_quadrilateral_on_the_same_plane(
        &nodes.get(&node_1_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
        &nodes.get(&node_2_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
        &nodes.get(&node_3_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
        &nodes.get(&node_4_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
        abs_tol,
    )
    {
        return Err(PlateElementDataError::<V>::NodesNotOnPlane(number).compose_error_message());
    }
    if convex_hull_on_four_points_on_plane(
            &[node_1_number, node_2_number, node_3_number, node_4_number],
            &[
                &nodes.get(&node_1_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
                &nodes.get(&node_2_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
                &nodes.get(&node_3_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
                &nodes.get(&node_4_number).ok_or(format!("Node {node_1_number} is absent!"))?.get_coordinates(),
            ],
            rel_tol,
            abs_tol,
        )?
        .len() != 4
    {
        return Err(PlateElementDataError::<V>::NotConvex(number).compose_error_message());
    }

    Ok(())
}


// fn power_func_x<V>(a: V, x: V, n: i32) -> V
//     where V: FloatTrait<Output = V>
// {
//     (0..n).fold(a, |acc, _| acc * x)
// }


// fn derivative_x<V>(f: fn(V, V, i32) -> V, a: V, x: V, n: i32) -> V
//     where V: FloatTrait<Output = V>
// {
//     let mut converted_n = V::from(0f32);
//     (0..n).for_each(|_| converted_n += V::from(1f32));
//     f(a * converted_n, x, n - 1)
// }


// fn dx_dr<V>(x_1: V, x_2: V, r: V) -> V
//     where V: FloatTrait<Output = V>
// {
//     derivative_x(power_func_x, x_1 * V::from(0.5f32), V::from(0f32), 0) -
//     derivative_x(power_func_x, x_1 * V::from(0.5f32), r, 1) +
//     derivative_x(power_func_x, x_2 * V::from(0.5f32), V::from(0f32), 0) +
//     derivative_x(power_func_x, x_2 * V::from(0.5f32), r, 1)
// }


// fn jacobian_at_r<V>(node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, Node<V>>) -> Result<V, String>
//     where V: FloatTrait<Output = V>
// {
//     let beam_element_vector = find_beam_element_vector(node_1_number, node_2_number, nodes)?;
//     let beam_element_length = beam_element_vector.norm()?;
//     let x_1 = V::from(-1f32) * beam_element_length / V::from(2f32);
//     let x_2 = beam_element_length / V::from(2f32);
//     Ok(dx_dr(x_1, x_2, r))
// }


// fn inverse_jacobian_at_r<V>(
//     node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, Node<V>>,
// ) 
//     -> Result<V, String>
//     where V: FloatTrait<Output = V>
// {
//     Ok(V::from(1f32) / jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
// }


// fn determinant_of_jacobian_at_r<V>(
//     node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, Node<V>>,
// ) 
//     -> Result<V, String>
//     where V: FloatTrait<Output = V>
// {
//     jacobian_at_r(node_1_number, node_2_number, r, nodes)
// }


// pub fn h1_r<V>(r: V) -> V
//     where V: FloatTrait<Output = V>
// {
//     V::from(0.5f32) * (V::from(1f32) - r)
// }


// pub fn h2_r<V>(r: V) -> V
//     where V: FloatTrait<Output = V>
// {
//     V::from(0.5f32) * (V::from(1f32) + r)
// }


// fn dh1_dr<V>(r: V) -> V
//     where V: FloatTrait<Output = V>
// {
//     derivative_x(power_func_x, V::from(0.5f32), V::from(0f32), 0) -
//     derivative_x(power_func_x, V::from(0.5f32), r, 1)
// }


// fn dh2_dr<V>(r: V) -> V
//     where V: FloatTrait<Output = V>
// {
//     derivative_x(power_func_x, V::from(0.5f32), V::from(0f32), 0) +
//     derivative_x(power_func_x, V::from(0.5f32), r, 1)
// }


// pub fn strain_displacement_matrix_u_at_r<V>(
//     node_1_number: u32, 
//     node_2_number: u32, 
//     r: V, 
//     nodes: &HashMap<u32, Node<V>>
// )
//     -> Result<Matrix<V>, String>
//     where V: FloatTrait<Output = V>
// {
//     let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;
//     Ok
//     (
//         Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[
//                 dh1_dr(r), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
//                 dh2_dr(r), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
//             ],
//         )
//         .multiply_by_scalar(inverse_jacobian)
//     )
// }


// fn strain_displacement_matrix_v_at_r<V>(
//     node_1_number: u32,
//     node_2_number: u32, 
//     r: V,
//     nodes: &HashMap<u32, Node<V>>,
// ) 
//     -> Result<Matrix<V>, String>
//     where V: FloatTrait<Output = V>
// {
//     let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;
    
//     let lhs_matrix = Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[
//                 V::from(0f32), dh1_dr(r), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
//                 V::from(0f32), dh2_dr(r), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
//             ],
//         )
//         .multiply_by_scalar(inverse_jacobian);

//     let rhs_matrix = Matrix::create(
//         1,
//         BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//         &[
//             V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), h1_r(r),
//             V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), h2_r(r),
//         ],
//     );

//     let matrix = lhs_matrix.subtract(&rhs_matrix)?;

//     Ok(matrix)
// }


// fn strain_displacement_matrix_w_at_r<V>(
//     node_1_number: u32, 
//     node_2_number: u32, 
//     r: V,
//     nodes: &HashMap<u32, Node<V>>
// ) 
//     -> Result<Matrix<V>, String>
//     where V: FloatTrait<Output = V>
// {
//     let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;

//     let lhs_matrix = Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[
//                 V::from(0f32), V::from(0f32), dh1_dr(r), V::from(0f32), V::from(0f32), V::from(0f32),
//                 V::from(0f32), V::from(0f32), dh2_dr(r), V::from(0f32), V::from(0f32), V::from(0f32),
//             ],
//         )
//         .multiply_by_scalar(inverse_jacobian);

//     let rhs_matrix = Matrix::create(
//         1,
//         BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//         &[
//             V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), h1_r(r), V::from(0f32),
//             V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), h2_r(r), V::from(0f32),
//         ],
//     );

//     let matrix = lhs_matrix.subtract(&rhs_matrix)?;

//     Ok(matrix)
// }


// fn strain_displacement_matrix_thu_at_r<V>(
//     node_1_number: u32, 
//     node_2_number: u32, 
//     r: V,
//     nodes: &HashMap<u32, Node<V>>,
// ) 
//     -> Result<Matrix<V>, String>
//     where V: FloatTrait<Output = V>
// {
//     let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;
//     Ok
//     (
//         Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[
//                 V::from(0f32), V::from(0f32), V::from(0f32), dh1_dr(r), V::from(0f32), V::from(0f32),
//                 V::from(0f32), V::from(0f32), V::from(0f32), dh2_dr(r), V::from(0f32), V::from(0f32),
//             ],
//         )
//         .multiply_by_scalar(inverse_jacobian)
//     )
// }


// fn strain_displacement_matrix_thv_at_r<V>(
//     node_1_number: u32, 
//     node_2_number: u32, 
//     r: V, 
//     nodes: &HashMap<u32, Node<V>>
// ) 
//     -> Result<Matrix<V>, String>
//     where V: FloatTrait<Output = V>
// {
//     let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;
//     Ok
//     (
//         Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[
//                 V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), dh1_dr(r), V::from(0f32),
//                 V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), dh2_dr(r), V::from(0f32),
//             ],
//         )
//         .multiply_by_scalar(inverse_jacobian)
//     )
// }


// pub fn strain_displacement_matrix_thw_at_r<V>(
//     node_1_number: u32, 
//     node_2_number: u32, 
//     r: V,
//     nodes: &HashMap<u32, Node<V>>,
// ) 
//     -> Result<Matrix<V>, String>
//     where V: FloatTrait<Output = V>
// {
//     let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;

//     Ok
//     (
//         Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[
//                 V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), dh1_dr(r),
//                 V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), dh2_dr(r),
//             ],
//         )
//         .multiply_by_scalar(inverse_jacobian)
//     )
// }


// fn local_stiffness_matrix_at_ip<V>(
//     node_1_number: u32,
//     node_2_number: u32,
//     young_modulus: V,
//     poisson_ratio: V,
//     area: V,
//     i11_p: V,
//     i22_p: V,
//     it: V,
//     shear_factor: V,
//     r: V,
//     alpha: V,
//     nodes: &HashMap<u32, Node<V>>,
// ) 
//     -> Result<SquareMatrix<V>, String>
//     where V: FloatTrait<Output = V>
// {
//     let b_u_at_r = strain_displacement_matrix_u_at_r(node_1_number, node_2_number, r, nodes)?;
//     let b_u_t_at_r = b_u_at_r.transpose();
//     let c_u_at_r = area * young_modulus;
//     let k_u_at_ip = b_u_t_at_r
//         .multiply_by_scalar(c_u_at_r)
//         .multiply(&b_u_at_r)?
//         .multiply_by_scalar(determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
//         .multiply_by_scalar(alpha);

//     let shear_modulus = young_modulus / (V::from(2f32) * (V::from(1f32) + poisson_ratio));

//     let b_v_at_r = strain_displacement_matrix_v_at_r(node_1_number, node_2_number, r, nodes)?;
//     let b_v_t_at_r = b_v_at_r.transpose();
//     let c_v_at_r = shear_modulus * area * shear_factor;
//     let k_v_at_ip = b_v_t_at_r
//         .multiply_by_scalar(c_v_at_r)
//         .multiply(&b_v_at_r)?
//         .multiply_by_scalar(determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
//         .multiply_by_scalar(alpha);

//     let b_w_at_r = strain_displacement_matrix_w_at_r(node_1_number, node_2_number, r, nodes)?;
//     let b_w_t_at_r = b_w_at_r.transpose();
//     let c_w_at_r = shear_modulus * area * shear_factor;
//     let k_w_at_ip = b_w_t_at_r
//         .multiply_by_scalar(c_w_at_r)
//         .multiply(&b_w_at_r)?
//         .multiply_by_scalar(determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
//         .multiply_by_scalar(alpha);

//     let b_thu_at_r = strain_displacement_matrix_thu_at_r(node_1_number, node_2_number, r, nodes)?;
//     let b_thu_t_at_r = b_thu_at_r.transpose();
//     let c_thu_at_r = shear_modulus * it;
//     let k_thu_at_ip = b_thu_t_at_r
//         .multiply_by_scalar(c_thu_at_r)
//         .multiply(&b_thu_at_r)?
//         .multiply_by_scalar(determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
//         .multiply_by_scalar(alpha);

//     let b_thv_at_r = strain_displacement_matrix_thv_at_r(node_1_number, node_2_number, r, nodes)?;
//     let b_thv_t_at_r = b_thv_at_r.transpose();
//     let c_thv_at_r = young_modulus * i22_p;
//     let k_thv_at_ip = b_thv_t_at_r
//         .multiply_by_scalar(c_thv_at_r)
//         .multiply(&b_thv_at_r)?
//         .multiply_by_scalar(determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
//         .multiply_by_scalar(alpha);

//     let b_thw_at_r = strain_displacement_matrix_thw_at_r(node_1_number, node_2_number, r, nodes)?;
//     let b_thw_t_at_r = b_thw_at_r.transpose();
//     let c_thw_at_r = young_modulus * i11_p;
//     let k_thw_at_ip = b_thw_t_at_r
//         .multiply_by_scalar(c_thw_at_r)
//         .multiply(&b_thw_at_r)?
//         .multiply_by_scalar(determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
//         .multiply_by_scalar(alpha);
    
//     Ok(
//         k_u_at_ip
//             .add(&k_v_at_ip)?
//             .add(&k_w_at_ip)?
//             .add(&k_thu_at_ip)?
//             .add(&k_thv_at_ip)?
//             .add(&k_thw_at_ip)?
//             .try_into_square_matrix()?
//     )
// }


// fn compose_local_stiffness_matrix<V>(
//     integration_points: &[(V, V)], 
//     node_1_number: u32,
//     node_2_number: u32,
//     young_modulus: V,
//     poisson_ratio: V,
//     area: V,
//     i11_p: V,
//     i22_p: V,
//     it: V,
//     shear_factor: V,
//     nodes: &HashMap<u32, Node<V>>,
// ) 
//     -> Result<SquareMatrix<V>, String>
//     where V: FloatTrait<Output = V>
// {
//     let mut local_stiffness_matrix = SquareMatrix::create(
//         BEAM_NODES_NUMBER * BEAM_NODE_DOF, 
//         &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
//     );

//     for (r, alpha) in integration_points
//     {
//         let local_stiffness_matrix_at_ip = local_stiffness_matrix_at_ip(
//             node_1_number,
//             node_2_number,
//             young_modulus,
//             poisson_ratio,
//             area,
//             i11_p,
//             i22_p,
//             it,
//             shear_factor,
//             *r,
//             *alpha,
//             nodes,
//         )?;
//         local_stiffness_matrix = local_stiffness_matrix.add(&local_stiffness_matrix_at_ip)?;
//     }

//     Ok(local_stiffness_matrix)
// }


fn compose_rotation_matrix<V>(rotation_matrix_elements: &[V; 9]) -> SquareMatrix<V>
    where V: FloatTrait<>
{
    let [q_11, q_12, q_13, q_21, q_22, q_23, q_31, q_32, q_33] = 
        rotation_matrix_elements;

    SquareMatrix::create(
        PLATE_NODES_NUMBER * PLATE_NODE_DOF,
        &[
            [*q_11, *q_12, *q_13, V::from(0f32), V::from(0f32), V::from(0f32)],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF],

            [*q_21, *q_22, *q_23, V::from(0f32), V::from(0f32), V::from(0f32)], 
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],

            [*q_31, *q_32, *q_33, V::from(0f32), V::from(0f32), V::from(0f32)], 
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32), V::from(0f32), V::from(0f32), *q_11, *q_12, *q_13],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32), V::from(0f32), V::from(0f32), *q_21, *q_22, *q_23],
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32), V::from(0f32), V::from(0f32), *q_31, *q_32, *q_33],
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],


            [V::from(0f32); PLATE_NODE_DOF], 
            [*q_11, *q_12, *q_13, V::from(0f32), V::from(0f32), V::from(0f32)],
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF],
            [*q_21, *q_22, *q_23, V::from(0f32), V::from(0f32), V::from(0f32)], 
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF],
            [*q_31, *q_32, *q_33, V::from(0f32), V::from(0f32), V::from(0f32)], 
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32), V::from(0f32), V::from(0f32), *q_11, *q_12, *q_13],
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32), V::from(0f32), V::from(0f32), *q_21, *q_22, *q_23],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32), V::from(0f32), V::from(0f32), *q_31, *q_32, *q_33],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],


            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [*q_11, *q_12, *q_13, V::from(0f32), V::from(0f32), V::from(0f32)],
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [*q_21, *q_22, *q_23, V::from(0f32), V::from(0f32), V::from(0f32)], 
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [*q_31, *q_32, *q_33, V::from(0f32), V::from(0f32), V::from(0f32)], 
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32), V::from(0f32), V::from(0f32), *q_11, *q_12, *q_13],
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32), V::from(0f32), V::from(0f32), *q_21, *q_22, *q_23],
            [V::from(0f32); PLATE_NODE_DOF],

            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32), V::from(0f32), V::from(0f32), *q_31, *q_32, *q_33],
            [V::from(0f32); PLATE_NODE_DOF],


            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [*q_11, *q_12, *q_13, V::from(0f32), V::from(0f32), V::from(0f32)],

            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [*q_21, *q_22, *q_23, V::from(0f32), V::from(0f32), V::from(0f32)], 

            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [*q_31, *q_32, *q_33, V::from(0f32), V::from(0f32), V::from(0f32)], 

            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32), V::from(0f32), V::from(0f32), *q_11, *q_12, *q_13],
            
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32), V::from(0f32), V::from(0f32), *q_21, *q_22, *q_23],
            
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32); PLATE_NODE_DOF], 
            [V::from(0f32); PLATE_NODE_DOF],
            [V::from(0f32), V::from(0f32), V::from(0f32), *q_31, *q_32, *q_33],
        ]
        .concat(),
    )
}


pub struct Plate<V>
{
    node_1_number: u32,
    node_2_number: u32,
    node_3_number: u32,
    node_4_number: u32,
    young_modulus: V,
    poisson_ratio: V,
    thickness: V,
    shear_factor: V,
    rotation_matrix_elements: [V; 9],
    integration_points: [(V, V, V, V); 4],
}


impl<V> Plate<V>
    where V: FloatTrait<Output = V>
{
    pub fn create(
        number: u32,
        node_1_number: u32,
        node_2_number: u32,
        node_3_number: u32,
        node_4_number: u32,
        young_modulus: V,
        poisson_ratio: V,
        thickness: V,
        shear_factor: V,
        nodes: &HashMap<u32, Node<V>>,
        rel_tol: V,
        abs_tol: V,
    )
        -> Result<Self, String>
    {
        check_plate_properties(
            number,
            young_modulus,
            poisson_ratio,
            thickness,
            shear_factor,
            node_1_number,
            node_2_number,
            node_3_number,
            node_4_number,
            nodes,
            rel_tol,
            abs_tol,
        )?;

        let rotation_matrix_elements = find_rotation_matrix_elements_of_quadrilateral(

            &nodes.get(&node_2_number)
                .ok_or(format!("Node {node_1_number} is absent!"))?
                .get_coordinates(),
            &nodes.get(&node_3_number)
                .ok_or(format!("Node {node_1_number} is absent!"))?
                .get_coordinates(),
            &nodes.get(&node_4_number)
                .ok_or(format!("Node {node_1_number} is absent!"))?
                .get_coordinates(),
            rel_tol,
            abs_tol,
        )?;
        let integration_points = [
            (
                V::from(1f32 / 3f32).my_sqrt() * V::from(1f32), V::from(1f32 / 3f32).my_sqrt() * V::from(1f32),
                V::from(1f32), V::from(1f32),
            ),
            (
                V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32), V::from(1f32 / 3f32).my_sqrt() * V::from(1f32),
                V::from(1f32), V::from(1f32),
            ),
            (
                V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32), V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32),
                V::from(1f32), V::from(1f32),
            ),
            (
                V::from(1f32 / 3f32).my_sqrt() * V::from(1f32), V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32),
                V::from(1f32), V::from(1f32),
            ),
        ];

        Ok(
            Plate 
            { 
                node_1_number, node_2_number, node_3_number, node_4_number, young_modulus, poisson_ratio, thickness, 
                shear_factor, rotation_matrix_elements, integration_points,
            }
        )
    }


    fn is_node_belongs_to_element(&self, node_number: u32) -> bool
    {
        self.node_1_number == node_number || self.node_2_number == node_number ||
        self.node_3_number == node_number || self.node_4_number == node_number
    }


    pub fn is_nodes_numbers_same(&self, nodes_numbers: &[u32; 4]) -> bool
    {
        nodes_numbers.iter().all(|node_number| self.is_node_belongs_to_element(*node_number))
    }


    pub fn extract_rotation_matrix(&self) -> SquareMatrix<V>
    {
        compose_rotation_matrix(&self.rotation_matrix_elements)
    }


    // pub fn extract_local_stiffness_matrix(&self, nodes: &HashMap<u32, Node<V>>) -> Result<SquareMatrix<V>, String>
    // {
    //     compose_local_stiffness_matrix(
    //         &self.integration_points,
    //         self.node_1_number,
    //         self.node_2_number,
    //         self.young_modulus,
    //         self.poisson_ratio,
    //         self.area,
    //         self.i11_p,
    //         self.i22_p,
    //         self.it,
    //         self.shear_factor,
    //         nodes,
    //     )
    // }


//     pub fn convert_uniformly_distributed_line_force_to_nodal_forces(
//         &self,
//         uniformly_distributed_line_force_value: V,
//         nodes: &HashMap<u32, Node<V>>,
//     ) 
//         -> Result<Vector<V>, String>
//     {
//         let distributed_force_matrix = Matrix::create(
//             1, 
//             1,
//             &[uniformly_distributed_line_force_value],
//         );

//         let mut nodal_forces = Vector::create(&[V::from(0f32); 2]);
//         for (r, alpha) in self.integration_points.iter()
//         {
//             let determinant_of_jacobian_at_r = determinant_of_jacobian_at_r(
//                 self.node_1_number, self.node_2_number, *r, nodes,
//             )?;

//             let displacement_interpolation_matrix = Vector::create(&[h1_r(*r), h2_r(*r)]);

//             let matrix = displacement_interpolation_matrix
//                 .multiply(&distributed_force_matrix)?
//                 .multiply_by_scalar(determinant_of_jacobian_at_r)
//                 .multiply_by_scalar(*alpha);
//             nodal_forces = nodal_forces.add(&matrix)?;
//         }

//         Ok(nodal_forces)
//     }


//     pub fn get_nodes_numbers(&self) -> [u32; 2]
//     {
//         [self.node_1_number, self.node_2_number]
//     }


//     pub fn extract_element_analysis_result(
//         &self, nodes: &HashMap<u32, Node<V>>, displacements: &Vector<V>,
//     )
//         -> Result<Vec<(ElementForceComponent, V)>, String>
//     {
//         let node_1_index = nodes.get(&self.node_1_number)
//             .ok_or(format!("Node: {} is absent!", self.node_1_number))?
//             .get_index();
//         let node_2_index = nodes.get(&self.node_2_number)
//             .ok_or(format!("Node: {} is absent!", self.node_2_number))?
//             .get_index();
//         let mut global_displacements = Vector::create(
//             &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
//         );

//         for i in 0..BEAM_NODE_DOF
//         {
//             *global_displacements.get_mut_element_value(&Position(i, 0))? = 
//                 *displacements.get_element_value(&Position(node_1_index * NODE_DOF + i, 0))?;
//         }

//         for i in 0..BEAM_NODE_DOF
//         {
//             *global_displacements.get_mut_element_value(&Position(i + BEAM_NODE_DOF, 0))? = 
//                 *displacements.get_element_value(&Position(node_2_index * NODE_DOF + i, 0))?;
//         }

//         let rotation_matrix = compose_rotation_matrix(&self.rotation_matrix_elements);
//         let local_displacements = rotation_matrix.multiply(&global_displacements)?;

//         let mut strain_displacement_matrix_u = Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
//         );
//         for (r, _alpha) in self.integration_points.iter()
//         {
//             let strain_displacement_matrix_u_at_r = strain_displacement_matrix_u_at_r(
//                 self.node_1_number, self.node_2_number, *r, nodes,
//             )?;
//             strain_displacement_matrix_u = strain_displacement_matrix_u.add(&strain_displacement_matrix_u_at_r)?;
//         }
//         let element_strains_u = strain_displacement_matrix_u.multiply(&local_displacements)?;
//         let element_forces_u = element_strains_u.multiply_by_scalar(
//             self.young_modulus * self.area / V::from(self.integration_points.len() as f32),
//         );

//         let shear_modulus = self.young_modulus / (V::from(2f32) * (V::from(1f32) + self.poisson_ratio));

//         let mut strain_displacement_matrix_v = Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
//         );
//         for (r, _alpha) in self.integration_points.iter()
//         {
//             let strain_displacement_matrix_v_at_r = strain_displacement_matrix_v_at_r(
//                 self.node_1_number, self.node_2_number, *r, nodes,
//             )?;
//             strain_displacement_matrix_v = strain_displacement_matrix_v.add(&strain_displacement_matrix_v_at_r)?;
//         }
//         let element_strains_v = strain_displacement_matrix_v.multiply(&local_displacements)?;
//         let element_forces_v = element_strains_v.multiply_by_scalar(
//             shear_modulus * self.area  * self.shear_factor / V::from(self.integration_points.len() as f32),
//         );
//         let force_s_value = *element_forces_v.get_element_value(&Position(0, 0))?;

//         let mut strain_displacement_matrix_w = Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
//         );
//         for (r, _alpha) in self.integration_points.iter()
//         {
//             let strain_displacement_matrix_w_at_r = strain_displacement_matrix_w_at_r(
//                 self.node_1_number, self.node_2_number, *r, nodes,
//             )?;
//             strain_displacement_matrix_w = strain_displacement_matrix_w.add(&strain_displacement_matrix_w_at_r)?;
//         }
//         let element_strains_w = strain_displacement_matrix_w.multiply(&local_displacements)?;
//         let element_forces_w = element_strains_w.multiply_by_scalar(
//             shear_modulus * self.area  * self.shear_factor / V::from(self.integration_points.len() as f32),
//         );
//         let force_t_value = *element_forces_w.get_element_value(&Position(0, 0))?;

//         let mut strain_displacement_matrix_thu = Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
//         );
//         for (r, _alpha) in self.integration_points.iter()
//         {
//             let strain_displacement_matrix_thu_at_r = strain_displacement_matrix_thu_at_r(
//                 self.node_1_number, self.node_2_number, *r, nodes,
//             )?;
//             strain_displacement_matrix_thu = strain_displacement_matrix_thu.add(&strain_displacement_matrix_thu_at_r)?;
//         }
//         let element_strains_thu = strain_displacement_matrix_thu.multiply(&local_displacements)?;
//         let element_forces_thu = element_strains_thu.multiply_by_scalar(
//             shear_modulus * self.it / V::from(self.integration_points.len() as f32),
//         );

//         let beam_element_vector = find_beam_element_vector(self.node_1_number, self.node_2_number, nodes)?;
//         let beam_element_length = beam_element_vector.norm()?;

//         let mut strain_displacement_matrix_thv = Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
//         );
//         for (r, _alpha) in self.integration_points.iter()
//         {
//             let strain_displacement_matrix_thv_at_r = strain_displacement_matrix_thv_at_r(
//                 self.node_1_number, self.node_2_number, *r, nodes,
//             )?;
//             strain_displacement_matrix_thv = strain_displacement_matrix_thv.add(&strain_displacement_matrix_thv_at_r)?;
//         }
//         let element_strains_thv = strain_displacement_matrix_thv.multiply(&local_displacements)?;
//         let element_forces_thv = element_strains_thv.multiply_by_scalar(
//             self.young_modulus * self.i22_p / V::from(self.integration_points.len() as f32),
//         );
//         let moment_s_average = *element_forces_thv.get_element_value(&Position(0, 0))?;
//         let moment_s_at_node_1 = moment_s_average + beam_element_length * force_t_value / V::from(2f32);
//         let moment_s_at_node_2 = moment_s_average - beam_element_length * force_t_value / V::from(2f32);

//         let mut strain_displacement_matrix_thw = Matrix::create(
//             1,
//             BEAM_NODES_NUMBER * BEAM_NODE_DOF,
//             &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
//         );
//         for (r, _alpha) in self.integration_points.iter()
//         {
//             let strain_displacement_matrix_thw_at_r = strain_displacement_matrix_thw_at_r(
//                 self.node_1_number, self.node_2_number, *r, nodes,
//             )?;
//             strain_displacement_matrix_thw = strain_displacement_matrix_thw.add(&strain_displacement_matrix_thw_at_r)?;
//         }
//         let element_strains_thw = strain_displacement_matrix_thw.multiply(&local_displacements)?;
//         let element_forces_thw = element_strains_thw.multiply_by_scalar(
//             self.young_modulus * self.i11_p / V::from(self.integration_points.len() as f32),
//         );
//         let moment_t_average = *element_forces_thw.get_element_value(&Position(0, 0))?;
//         let moment_t_at_node_1 = moment_t_average + beam_element_length * force_s_value / V::from(2f32);
//         let moment_t_at_node_2 = moment_t_average - beam_element_length * force_s_value / V::from(2f32);

//         let element_analysis_data = vec![
//             (
//                 ElementForceComponent::ForceR,
//                 *element_forces_u.get_element_value(&Position(0, 0))?,
//             ),
//             (
//                 ElementForceComponent::ForceS,
//                 force_s_value,
//             ),
//             (
//                 ElementForceComponent::ForceT,
//                 force_t_value,
//             ),
//             (
//                 ElementForceComponent::MomentR,
//                 *element_forces_thu.get_element_value(&Position(0, 0))?,
//             ),
//             (
//                 ElementForceComponent::MomentS,
//                 moment_s_at_node_1,
//             ),
//             (
//                 ElementForceComponent::MomentS,
//                 moment_s_average,
//             ),
//             (
//                 ElementForceComponent::MomentS,
//                 moment_s_at_node_2,
//             ),
//             (
//                 ElementForceComponent::MomentT,
//                 moment_t_at_node_1,
//             ),
//             (
//                 ElementForceComponent::MomentT,
//                 moment_t_average,
//             ),
//             (
//                 ElementForceComponent::MomentT,
//                 moment_t_at_node_2,
//             ),
//         ];

//         Ok(element_analysis_data)
//     }
}
