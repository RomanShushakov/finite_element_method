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
    find_rotation_matrix_elements_of_quadrilateral, convex_hull_on_four_points_on_plane, dh_dx_dh_dy,
    extract_transformed_directions_of_nodes, determinant_of_jacobian_at_r_s,
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


fn strain_displacement_matrix_mem_at_r_s<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    node_3_number: u32, 
    node_4_number: u32,
    r: V, 
    s: V, 
    nodes: &HashMap<u32, Node<V>>, 
    rotation_matrix_elements: &[V; 9],
    rel_tol: V,
) 
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let dh_dx_dh_dy_matrix = dh_dx_dh_dy(
        node_1_number, node_2_number, node_3_number, node_4_number, r, s, nodes, rotation_matrix_elements, rel_tol,
    )?;

    let dh1_dx = *dh_dx_dh_dy_matrix.get_element_value(&Position(0, 0))?;
    let dh2_dx = *dh_dx_dh_dy_matrix.get_element_value(&Position(0, 1))?;
    let dh3_dx = *dh_dx_dh_dy_matrix.get_element_value(&Position(0, 2))?;
    let dh4_dx = *dh_dx_dh_dy_matrix.get_element_value(&Position(0, 3))?;

    let dh1_dy = *dh_dx_dh_dy_matrix.get_element_value(&Position(1, 0))?;
    let dh2_dy = *dh_dx_dh_dy_matrix.get_element_value(&Position(1, 1))?;
    let dh3_dy = *dh_dx_dh_dy_matrix.get_element_value(&Position(1, 2))?;
    let dh4_dy = *dh_dx_dh_dy_matrix.get_element_value(&Position(1, 3))?;

    Ok(
        Matrix::create(
            3, 
            PLATE_NODES_NUMBER * PLATE_NODE_DOF, 
            &[
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
            ],
        )
    )
}


fn strain_displacement_matrix_bend_at_r_s<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    node_3_number: u32, 
    node_4_number: u32,
    r: V, 
    s: V, 
    nodes: &HashMap<u32, Node<V>>, 
    rotation_matrix_elements: &[V; 9],
    rel_tol: V,
) 
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let dh_dx_dh_dy_matrix = dh_dx_dh_dy(
        node_1_number, node_2_number, node_3_number, node_4_number, r, s, nodes, rotation_matrix_elements, rel_tol,
    )?;

    let dh1_dx = *dh_dx_dh_dy_matrix.get_element_value(&Position(0, 0))?;
    let dh2_dx = *dh_dx_dh_dy_matrix.get_element_value(&Position(0, 1))?;
    let dh3_dx = *dh_dx_dh_dy_matrix.get_element_value(&Position(0, 2))?;
    let dh4_dx = *dh_dx_dh_dy_matrix.get_element_value(&Position(0, 3))?;

    let dh1_dy = *dh_dx_dh_dy_matrix.get_element_value(&Position(1, 0))?;
    let dh2_dy = *dh_dx_dh_dy_matrix.get_element_value(&Position(1, 1))?;
    let dh3_dy = *dh_dx_dh_dy_matrix.get_element_value(&Position(1, 2))?;
    let dh4_dy = *dh_dx_dh_dy_matrix.get_element_value(&Position(1, 3))?;

    Ok(
        Matrix::create(
            3, 
            PLATE_NODES_NUMBER * PLATE_NODE_DOF, 
            &[
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
            ],
        )
    )
}


fn strain_displacement_matrix_shear_at_r_s<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    node_3_number: u32, 
    node_4_number: u32,
    r: V, 
    s: V, 
    nodes: &HashMap<u32, Node<V>>, 
    rotation_matrix_elements: &[V; 9],
    rel_tol: V,
) 
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let transformed_directions_of_nodes = extract_transformed_directions_of_nodes(
        node_1_number, node_2_number, node_3_number, node_4_number, nodes, rotation_matrix_elements,
    )?;

    let x_1 = transformed_directions_of_nodes[0][0];
    let y_1 = transformed_directions_of_nodes[0][1];
    let x_2 = transformed_directions_of_nodes[1][0];
    let y_2 = transformed_directions_of_nodes[1][1];
    let x_3 = V::from(0f32);
    let y_3 = V::from(0f32);
    let x_4 = transformed_directions_of_nodes[2][0];
    let y_4 = transformed_directions_of_nodes[2][1];

    let a_x = x_1 - x_2 - x_3 + x_4;
    let b_x = x_1 - x_2 + x_3 - x_4;
    let c_x = x_1 + x_2 - x_3 - x_4;
    let a_y = y_1 - y_2 - y_3 + y_4;
    let b_y = y_1 - y_2 + y_3 - y_4;
    let c_y = y_1 + y_2 - y_3 - y_4;

    let determinant_of_jacobian = determinant_of_jacobian_at_r_s(
        node_1_number, node_2_number, node_3_number, node_4_number, r, s, nodes, rotation_matrix_elements, rel_tol,
    )?;

    let gamma_rz_multiplier = (
        (c_x + r * b_x).my_powi(2) + (c_y + r * b_y).my_powi(2)).my_sqrt() / (V::from(8f32) * determinant_of_jacobian
    );
    
    let gamma_sz_multiplier = (
        (a_x + s * b_x).my_powi(2) + (a_y + s * b_y).my_powi(2)).my_sqrt() / (V::from(8f32) * determinant_of_jacobian
    );

    Ok(
        Matrix::create(
            2, 
            PLATE_NODES_NUMBER * PLATE_NODE_DOF, 
            &[
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
            ],
        )
    )
}


fn local_stiffness_matrix_at_ip<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    node_3_number: u32, 
    node_4_number: u32,
    young_modulus: V, 
    poisson_ratio: V, 
    thickness: V, 
    shear_factor: V,  
    r: V, 
    s: V,
    alpha: V,
    nodes: &HashMap<u32, Node<V>>, 
    rotation_matrix_elements: &[V; 9], 
    rel_tol: V,
) 
    -> Result<SquareMatrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let c_multiplier_mem = young_modulus * thickness / (V::from(1f32) - poisson_ratio.my_powi(2));
    let c_mem = Matrix::create(
            3, 
            3, 
            &[
                V::from(1f32), poisson_ratio, V::from(0f32),
                poisson_ratio, V::from(1f32), V::from(0f32),
                V::from(0f32), V::from(0f32), (V::from(1f32) - poisson_ratio) / V::from(2f32),
            ], 
        )
        .multiply_by_scalar(c_multiplier_mem);

    let b_mem_at_r_s = strain_displacement_matrix_mem_at_r_s(
        node_1_number, node_2_number, node_3_number, node_4_number, r, s, nodes, rotation_matrix_elements, rel_tol,
    )?;
    let b_mem_t_at_r_s = b_mem_at_r_s.transpose();

    let k_mem_at_ip = b_mem_t_at_r_s
        .multiply(&c_mem)?
        .multiply(&b_mem_at_r_s)?
        .multiply_by_scalar(
            determinant_of_jacobian_at_r_s(
                node_1_number, 
                node_2_number, 
                node_3_number, 
                node_4_number, 
                r, 
                s, 
                nodes, 
                rotation_matrix_elements, 
                rel_tol,
            )? * alpha,
        );

    let c_multiplier_bend = young_modulus * thickness.my_powi(3) / 
        (V::from(12f32) * (V::from(1f32) - poisson_ratio.my_powi(2)));
    let c_bend = Matrix::create(
            3, 
            3, 
            &[
                V::from(1f32), poisson_ratio, V::from(0f32),
                poisson_ratio, V::from(1f32), V::from(0f32),
                V::from(0f32), V::from(0f32), (V::from(1f32) - poisson_ratio) / V::from(2f32),
            ], 
        )
        .multiply_by_scalar(c_multiplier_bend);

    let b_bend_at_r_s = strain_displacement_matrix_bend_at_r_s(
        node_1_number, node_2_number, node_3_number, node_4_number, r, s, nodes, rotation_matrix_elements, rel_tol,
    )?;
    let b_bend_t_at_r_s = b_bend_at_r_s.transpose();

    let k_bend_at_ip = b_bend_t_at_r_s
        .multiply(&c_bend)?
        .multiply(&b_bend_at_r_s)?
        .multiply_by_scalar(
            determinant_of_jacobian_at_r_s(
                node_1_number, 
                node_2_number, 
                node_3_number, 
                node_4_number, 
                r, 
                s, 
                nodes, 
                rotation_matrix_elements,
                rel_tol
            )? * alpha,
        );

    let c_multiplier_shear = young_modulus * thickness * shear_factor / 
        (V::from(2f32) * (V::from(1f32) + poisson_ratio));
    let c_shear = Matrix::create(
            2, 
            2, 
            &[V::from(1f32), V::from(0f32), V::from(0f32), V::from(1f32)], 
        )
        .multiply_by_scalar(c_multiplier_shear * alpha);

    let b_shear_at_r_s = strain_displacement_matrix_shear_at_r_s(
        node_1_number, node_2_number, node_3_number, node_4_number, r, s, nodes, rotation_matrix_elements, rel_tol,
    )?;
    let b_shear_t_at_r_s = b_shear_at_r_s.transpose();

    let k_shear_at_ip = b_shear_t_at_r_s
        .multiply(&c_shear)?
        .multiply(&b_shear_at_r_s)?
        .multiply_by_scalar(
            determinant_of_jacobian_at_r_s(
                node_1_number, 
                node_2_number, 
                node_3_number, 
                node_4_number, 
                r, 
                s, 
                nodes, 
                rotation_matrix_elements,
                rel_tol,
            )? * alpha,
        );

    Ok(
        k_mem_at_ip
            .add(&k_bend_at_ip)?
            .add(&k_shear_at_ip)?
            .try_into_square_matrix()?
    )
}


fn compose_local_stiffness_matrix<V>(
    integration_points: &[(V, V, V, V); 4], 
    node_1_number: u32,
    node_2_number: u32,
    node_3_number: u32,
    node_4_number: u32,
    young_modulus: V,
    poisson_ratio: V,
    thickness: V,
    shear_factor: V,
    nodes: &HashMap<u32, Node<V>>,
    rotation_matrix_elements: &[V; 9],
    rel_tol: V,
) 
    -> Result<SquareMatrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let mut local_stiffness_matrix = SquareMatrix::create(
        PLATE_NODES_NUMBER * PLATE_NODE_DOF, 
        &[V::from(0f32); PLATE_NODES_NUMBER * PLATE_NODE_DOF],
    );

    for (r, s, alpha_r, alpha_s) in integration_points
    {
        let local_stiffness_matrix_at_ip = local_stiffness_matrix_at_ip(
            node_1_number,
            node_2_number,
            node_3_number,
            node_4_number,
            young_modulus,
            poisson_ratio,
            thickness,
            shear_factor,
            *r,
            *s,
            *alpha_r * *alpha_s,
            nodes,
            rotation_matrix_elements,
            rel_tol,
        )?;
        local_stiffness_matrix = local_stiffness_matrix.add(&local_stiffness_matrix_at_ip)?;
    }

    let krot6 = V::from(1f32);
    for i in 0..PLATE_NODES_NUMBER
    {
        *local_stiffness_matrix.get_mut_element_value(
            &Position(i * PLATE_NODE_DOF + PLATE_NODE_DOF - 1, i * PLATE_NODE_DOF + PLATE_NODE_DOF - 1),
        )? += krot6;
    }

    Ok(local_stiffness_matrix)
}


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


    pub fn extract_local_stiffness_matrix(
        &self, nodes: &HashMap<u32, Node<V>>, rel_tol: V,
    ) 
        -> Result<SquareMatrix<V>, String>
    {
        compose_local_stiffness_matrix(
            &self.integration_points,
            self.node_1_number,
            self.node_2_number,
            self.node_3_number,
            self.node_4_number,
            self.young_modulus,
            self.poisson_ratio,
            self.thickness,
            self.shear_factor,
            nodes,
            &self.rotation_matrix_elements,
            rel_tol,
        )
    }


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


    pub fn extract_element_analysis_result(
        &self, nodes: &HashMap<u32, Node<V>>, displacements: &Vector<V>, rel_tol: V,
    )
        -> Result<Vec<(ElementForceComponent, V)>, String>
    {
        let node_1_index = nodes.get(&self.node_1_number)
            .ok_or(format!("Node: {} is absent!", self.node_1_number))?
            .get_index();
        let node_2_index = nodes.get(&self.node_2_number)
            .ok_or(format!("Node: {} is absent!", self.node_2_number))?
            .get_index();
        let node_3_index = nodes.get(&self.node_3_number)
            .ok_or(format!("Node: {} is absent!", self.node_3_number))?
            .get_index();
        let node_4_index = nodes.get(&self.node_4_number)
            .ok_or(format!("Node: {} is absent!", self.node_4_number))?
            .get_index();
        let mut global_displacements = Vector::create(
            &[V::from(0f32); PLATE_NODES_NUMBER * PLATE_NODE_DOF],
        );

        for i in 0..PLATE_NODE_DOF
        {
            *global_displacements.get_mut_element_value(&Position(i, 0))? = 
                *displacements.get_element_value(&Position(node_1_index * NODE_DOF + i, 0))?;
        }

        for i in 0..PLATE_NODE_DOF
        {
            *global_displacements.get_mut_element_value(&Position(i + PLATE_NODE_DOF, 0))? = 
                *displacements.get_element_value(&Position(node_2_index * NODE_DOF + i, 0))?;
        }

        for i in 0..PLATE_NODE_DOF
        {
            *global_displacements.get_mut_element_value(&Position(i + PLATE_NODE_DOF * 2, 0))? = 
                *displacements.get_element_value(&Position(node_3_index * NODE_DOF + i, 0))?;
        }

        for i in 0..PLATE_NODE_DOF
        {
            *global_displacements.get_mut_element_value(&Position(i + PLATE_NODE_DOF * 3, 0))? = 
                *displacements.get_element_value(&Position(node_4_index * NODE_DOF + i, 0))?;
        }

        let rotation_matrix = compose_rotation_matrix(&self.rotation_matrix_elements);
        let local_displacements = rotation_matrix.multiply(&global_displacements)?;

        let local_nodes_coordinates = vec![
            (V::from(1f32), V::from(1f32)), 
            (V::from(-1f32), V::from(1f32)),
            (V::from(-1f32), V::from(-1f32)), 
            (V::from(1f32), V::from(-1f32))
        ];

        let c_multiplier_mem = self.young_modulus / (V::from(1f32) - self.poisson_ratio.my_powi(2));
        let c_mem = Matrix::create(
                3, 
                3, 
                &[
                    V::from(1f32), self.poisson_ratio, V::from(0f32),
                    self.poisson_ratio, V::from(1f32), V::from(0f32),
                    V::from(0f32), V::from(0f32), (V::from(1f32) - self.poisson_ratio) / V::from(2f32),
                ], 
            )
            .multiply_by_scalar(c_multiplier_mem);
        let mut strain_displacement_matrix_mem = Matrix::create(
            3,
            PLATE_NODES_NUMBER * PLATE_NODE_DOF,
            &[V::from(0f32); PLATE_NODES_NUMBER * PLATE_NODE_DOF * 3],
        );
        for (r, s) in local_nodes_coordinates.iter()
        {
            let strain_displacement_matrix_mem_at_r_s = strain_displacement_matrix_mem_at_r_s(
                self.node_1_number, 
                self.node_2_number,
                self.node_3_number,
                self.node_4_number,
                *r,
                *s,
                nodes,
                &self.rotation_matrix_elements,
                rel_tol,
            )?;
            strain_displacement_matrix_mem = strain_displacement_matrix_mem.add(&strain_displacement_matrix_mem_at_r_s)?;
        }
        let element_strains_mem = strain_displacement_matrix_mem.multiply(&local_displacements)?;
        let element_forces_mem = c_mem
            .multiply(&element_strains_mem)?
            .multiply_by_scalar(self.thickness / V::from(local_nodes_coordinates.len() as f32));

        let c_multiplier_bend = self.young_modulus * self.thickness / 
            (V::from(2f32) * (V::from(1f32) - self.poisson_ratio.my_powi(2)));
        let c_bend = Matrix::create(
                3, 
                3, 
                &[
                    V::from(1f32), self.poisson_ratio, V::from(0f32),
                    self.poisson_ratio, V::from(1f32), V::from(0f32),
                    V::from(0f32), V::from(0f32), (V::from(1f32) - self.poisson_ratio) / V::from(2f32),
                ], 
            )
            .multiply_by_scalar(c_multiplier_bend);
        let mut strain_displacement_matrix_bend = Matrix::create(
            3,
            PLATE_NODES_NUMBER * PLATE_NODE_DOF,
            &[V::from(0f32); PLATE_NODES_NUMBER * PLATE_NODE_DOF * 3],
        );
        for (r, s) in local_nodes_coordinates.iter()
        {
            let strain_displacement_matrix_bend_at_r_s = strain_displacement_matrix_bend_at_r_s(
                self.node_1_number, 
                self.node_2_number,
                self.node_3_number,
                self.node_4_number,
                *r,
                *s,
                nodes,
                &self.rotation_matrix_elements,
                rel_tol,
            )?;
            strain_displacement_matrix_bend = strain_displacement_matrix_bend.add(&strain_displacement_matrix_bend_at_r_s)?;
        }
        let element_strains_bend = strain_displacement_matrix_bend.multiply(&local_displacements)?;
        let element_forces_bend = c_bend
            .multiply(&element_strains_bend)?
            .multiply_by_scalar(self.thickness.my_powi(2) / V::from(24f32));

        let c_multiplier_shear = self.young_modulus / (V::from(2f32) * (V::from(1f32) + self.poisson_ratio));
        let c_shear = Matrix::create(
                2, 
                2, 
                &[
                    V::from(1f32), V::from(0f32),
                    V::from(0f32), V::from(1f32),
                ], 
            )
            .multiply_by_scalar(c_multiplier_shear);
        let mut strain_displacement_matrix_shear = Matrix::create(
            2,
            PLATE_NODES_NUMBER * PLATE_NODE_DOF,
            &[V::from(0f32); PLATE_NODES_NUMBER * PLATE_NODE_DOF * 2],
        );
        for (r, s) in local_nodes_coordinates.iter()
        {
            let strain_displacement_matrix_shear_at_r_s = strain_displacement_matrix_shear_at_r_s(
                self.node_1_number, 
                self.node_2_number,
                self.node_3_number,
                self.node_4_number,
                *r,
                *s,
                nodes,
                &self.rotation_matrix_elements,
                rel_tol,
            )?;
            strain_displacement_matrix_shear = strain_displacement_matrix_shear.add(&strain_displacement_matrix_shear_at_r_s)?;
        }
        let element_strains_shear = strain_displacement_matrix_shear.multiply(&local_displacements)?;
        let element_forces_shear = c_shear
            .multiply(&element_strains_shear)?
            .multiply_by_scalar(self.thickness * self.shear_factor / V::from(local_nodes_coordinates.len() as f32));

        let element_analysis_data = vec![
            (
                ElementForceComponent::MembraneForceR,
                *element_forces_mem.get_element_value(&Position(0, 0))?,
            ),
            (
                ElementForceComponent::MembraneForceS,
                *element_forces_mem.get_element_value(&Position(1, 0))?,
            ),
            (
                ElementForceComponent::MembraneForceRS,
                *element_forces_mem.get_element_value(&Position(2, 0))?,
            ),
            (
                ElementForceComponent::BendingMomentR,
                *element_forces_bend.get_element_value(&Position(1, 0))?,
            ),
            (
                ElementForceComponent::BendingMomentS,
                *element_forces_bend.get_element_value(&Position(0, 0))?,
            ),
            (
                ElementForceComponent::BendingMomentRS,
                *element_forces_bend.get_element_value(&Position(2, 0))?,
            ),
            (
                ElementForceComponent::ShearForceRT,
                *element_forces_shear.get_element_value(&Position(0, 0))?,
            ),
            (
                ElementForceComponent::ShearForceST,
                *element_forces_shear.get_element_value(&Position(1, 0))?,
            ),
        ];

        Ok(element_analysis_data)
    }
}
