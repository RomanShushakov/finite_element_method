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
use crate::fem::math_functions::compare_with_tolerance;
use crate::fem::bar_2n_element_functions::
{
    find_2n_element_vector, inverse_jacobian_at_r, determinant_of_jacobian_at_r, h1_r, h2_r, dh1_dr, dh2_dr,
};


const BEAM_NODES_NUMBER: usize = 2;
pub const BEAM_NODE_DOF: usize = 6;


enum BeamElementDataError<V>
{
    YoungModulus(V),
    PoissonRatio(V),
    Area(V),
    I11(V),
    I22(V),
    It(V),
    ShearFactor(V),
    ParallelLocalAxis(u32, V, V, V),
}


impl<V> BeamElementDataError<V>
    where V: Debug
{
    fn compose_error_message(&self) -> String
    {
        match self
        {
            Self::YoungModulus(value) => format!("Young's modulus {value:?} is less or equal to zero!"),
            Self::PoissonRatio(value) => format!("Poisson's ratio {value:?} is less or equal to zero!"),
            Self::Area(value) => format!("Area {value:?} is less or equal to zero!"),
            Self::I11(value) => format!("I11 {value:?} is less or equal to zero!"),
            Self::I22(value) => format!("I22 {value:?} is less or equal to zero!"),
            Self::It(value) => format!("It {value:?} is less or equal to zero!"),
            Self::ShearFactor(value) => format!("Shear factor {value:?} is less or equal to zero!"),
            Self::ParallelLocalAxis(number, x, y, z) => 
            {
                format!("Local axis 1 direction [{x:?}, {y:?}, {z:?}] parallel to element {number}!")
            },
        }
    }
}


fn check_beam_properties<V>(
    number: u32,
    young_modulus: V, 
    poisson_ratio: V, 
    area: V, 
    i11: V, 
    i22: V, 
    it: V, 
    shear_factor: V, 
    node_1_number: u32,
    node_2_number: u32,
    nodes: &HashMap<u32, Node<V>>,
    local_axis_1_direction: &[V; 3],
) 
    -> Result<(), String>
    where V: FloatTrait<Output = V>
{
    if young_modulus <= V::from(0f32)
    {
        return Err(BeamElementDataError::<V>::YoungModulus(young_modulus).compose_error_message());
    }
    if poisson_ratio <= V::from(0f32)
    {
        return Err(BeamElementDataError::<V>::PoissonRatio(young_modulus).compose_error_message());
    }
    if area <= V::from(0f32)
    {
        return Err(BeamElementDataError::<V>::Area(area).compose_error_message());
    }
    if i11 <= V::from(0f32)
    {
        return Err(BeamElementDataError::<V>::I11(i11).compose_error_message());
    }
    if i22 <= V::from(0f32)
    {
        return Err(BeamElementDataError::<V>::I22(i22).compose_error_message());
    }
    if it <= V::from(0f32)
    {
        return Err(BeamElementDataError::<V>::It(it).compose_error_message());
    }
    if shear_factor <= V::from(0f32)
    {
        return Err(BeamElementDataError::<V>::ShearFactor(shear_factor).compose_error_message());
    }

    let beam_element_vector = find_2n_element_vector(node_1_number, node_2_number, nodes)?;
    let local_axis_1_direction_vector = Vector3::create(local_axis_1_direction);
    let projection_of_beam_section_orientation = local_axis_1_direction_vector
        .projection_perpendicular_to_vector(&beam_element_vector);
    if projection_of_beam_section_orientation.norm()? == V::from(0f32)
    {
        return Err(BeamElementDataError::<V>::ParallelLocalAxis(
                number, local_axis_1_direction[0], local_axis_1_direction[1], local_axis_1_direction[2]
            )
            .compose_error_message()
        );
    }

    Ok(())
}


fn find_principal_moments_of_inertia<V>(i11: V, i22: V, i12: V) -> (V, V, V)
    where V: FloatTrait<Output = V>
{
    let mut angle = (V::from(2f32) * i12 / (i22 - i11)).my_atan() / V::from(2f32);

    let mut i11_p = i11 * (angle.my_cos()).my_powi(2) +
        i22 * (angle.my_sin()).my_powi(2) -
        i12 * (V::from(2f32) * angle).my_sin();

    let mut i22_p = i11 * (angle.my_sin()).my_powi(2) +
        i22 * (angle.my_cos()).my_powi(2) +
        i12 * (V::from(2f32) * angle).my_sin();

    let mut i = 1;
    while i11_p < i22_p
    {
        angle = ((V::from(2f32) * i12 / (i22 - i11)).my_atan() +
            V::from(PI) * V::from(i as f32)) / V::from(2f32);
        i11_p = i11 * (angle.my_cos()).my_powi(2) +
            i22 * (angle.my_sin()).my_powi(2) -
            i12 * (V::from(2f32) * angle).my_sin();
        i22_p = i11 * (angle.my_sin()).my_powi(2) +
            i22 * (angle.my_cos()).my_powi(2) +
            i12 * (V::from(2f32) * angle).my_sin();
        i += 1;
    }

    (i11_p, i22_p, angle)
}


fn find_rotation_matrix_elements<V>(
    node_1_number: u32,
    node_2_number: u32,
    local_axis_1_direction: &[V; 3],
    angle: V,
    nodes: &HashMap<u32, Node<V>>,
    rel_tol: V,
    abs_tol: V,
)
    -> Result<[V; 9], String>
    where V: FloatTrait<Output = V>
{
    let beam_element_vector = find_2n_element_vector(node_1_number, node_2_number, nodes)?;

    let beam_element_length = beam_element_vector.norm()?;
    let direction_vector = Vector3::create(
        &[beam_element_length, V::from(0f32), V::from(0f32)],
    );

    let interim_rotation_matrix = beam_element_vector
        .rotation_matrix_to_align_with_vector(&direction_vector, rel_tol, abs_tol)?;

    let local_axis_1_direction_vector = Vector3::create(local_axis_1_direction);
    let projection_of_beam_section_orientation = local_axis_1_direction_vector
        .projection_perpendicular_to_vector(&beam_element_vector);

    let transformed_projection_of_beam_section_orientation = interim_rotation_matrix
        .multiply(&projection_of_beam_section_orientation)?;

    let angle_between_beam_section_local_axis_1_direction_and_axis_t = Vector3::create(
        &[V::from(0f32), V::from(0f32), V::from(1f32)]
        )
        .cosine_angle_between_vectors(&Vector3::create(&[
            *transformed_projection_of_beam_section_orientation.get_element_value(&Position(0, 0))?,
            *transformed_projection_of_beam_section_orientation.get_element_value(&Position(1, 0))?,
            *transformed_projection_of_beam_section_orientation.get_element_value(&Position(2, 0))?,
        ]))
        .my_acos();

    let total_angle = angle + angle_between_beam_section_local_axis_1_direction_and_axis_t;
    let [x, y, z] = beam_element_vector.get_components();

    let c_x = compare_with_tolerance(x / beam_element_length, abs_tol);
    let c_y = compare_with_tolerance(y / beam_element_length, abs_tol);
    let c_z = compare_with_tolerance(z / beam_element_length, abs_tol);
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

    Ok([r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33])
}


fn strain_displacement_matrix_u_at_r<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    r: V, 
    nodes: &HashMap<u32, Node<V>>
)
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;
    Ok(
        Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[
                dh1_dr(r), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
                dh2_dr(r), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            ],
        )
        .multiply_by_scalar(inverse_jacobian)
    )
}


fn strain_displacement_matrix_v_at_r<V>(
    node_1_number: u32,
    node_2_number: u32, 
    r: V,
    nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;
    
    let lhs_matrix = Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[
                V::from(0f32), dh1_dr(r), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), dh2_dr(r), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32),
            ],
        )
        .multiply_by_scalar(inverse_jacobian);

    let rhs_matrix = Matrix::create(
        1,
        BEAM_NODES_NUMBER * BEAM_NODE_DOF,
        &[
            V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), h1_r(r),
            V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), h2_r(r),
        ],
    );

    let matrix = lhs_matrix.subtract(&rhs_matrix)?;

    Ok(matrix)
}


fn strain_displacement_matrix_w_at_r<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    r: V,
    nodes: &HashMap<u32, Node<V>>
) 
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;

    let lhs_matrix = Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[
                V::from(0f32), V::from(0f32), dh1_dr(r), V::from(0f32), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), dh2_dr(r), V::from(0f32), V::from(0f32), V::from(0f32),
            ],
        )
        .multiply_by_scalar(inverse_jacobian);

    let rhs_matrix = Matrix::create(
        1,
        BEAM_NODES_NUMBER * BEAM_NODE_DOF,
        &[
            V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), h1_r(r), V::from(0f32),
            V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), h2_r(r), V::from(0f32),
        ],
    );

    let matrix = lhs_matrix.subtract(&rhs_matrix)?;

    Ok(matrix)
}


fn strain_displacement_matrix_thu_at_r<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    r: V,
    nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;
    Ok(
        Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[
                V::from(0f32), V::from(0f32), V::from(0f32), dh1_dr(r), V::from(0f32), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32), dh2_dr(r), V::from(0f32), V::from(0f32),
            ],
        )
        .multiply_by_scalar(inverse_jacobian)
    )
}


fn strain_displacement_matrix_thv_at_r<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    r: V, 
    nodes: &HashMap<u32, Node<V>>
) 
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;
    Ok(
        Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[
                V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), dh1_dr(r), V::from(0f32),
                V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), dh2_dr(r), V::from(0f32),
            ],
        )
        .multiply_by_scalar(inverse_jacobian)
    )
}


fn strain_displacement_matrix_thw_at_r<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    r: V,
    nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<Matrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let inverse_jacobian = inverse_jacobian_at_r(node_1_number, node_2_number, r, nodes)?;

    Ok(
        Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[
                V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), dh1_dr(r),
                V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), V::from(0f32), dh2_dr(r),
            ],
        )
        .multiply_by_scalar(inverse_jacobian)
    )
}


fn local_stiffness_matrix_at_ip<V>(
    node_1_number: u32,
    node_2_number: u32,
    young_modulus: V,
    poisson_ratio: V,
    area: V,
    i11_p: V,
    i22_p: V,
    it: V,
    shear_factor: V,
    r: V,
    alpha: V,
    nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<SquareMatrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let b_u_at_r = strain_displacement_matrix_u_at_r(node_1_number, node_2_number, r, nodes)?;
    let b_u_t_at_r = b_u_at_r.transpose();
    let c_u = area * young_modulus;
    let k_u_at_ip = b_u_t_at_r
        .multiply(&b_u_at_r)?
        .multiply_by_scalar(c_u * determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)? * alpha);

    let shear_modulus = young_modulus / (V::from(2f32) * (V::from(1f32) + poisson_ratio));

    let b_v_at_r = strain_displacement_matrix_v_at_r(node_1_number, node_2_number, r, nodes)?;
    let b_v_t_at_r = b_v_at_r.transpose();
    let c_v = shear_modulus * area * shear_factor;
    let k_v_at_ip = b_v_t_at_r
        .multiply(&b_v_at_r)?
        .multiply_by_scalar(c_v * determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)? * alpha);

    let b_w_at_r = strain_displacement_matrix_w_at_r(node_1_number, node_2_number, r, nodes)?;
    let b_w_t_at_r = b_w_at_r.transpose();
    let c_w = shear_modulus * area * shear_factor;
    let k_w_at_ip = b_w_t_at_r
        .multiply(&b_w_at_r)?
        .multiply_by_scalar(c_w * determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)? * alpha);

    let b_thu_at_r = strain_displacement_matrix_thu_at_r(node_1_number, node_2_number, r, nodes)?;
    let b_thu_t_at_r = b_thu_at_r.transpose();
    let c_thu = shear_modulus * it;
    let k_thu_at_ip = b_thu_t_at_r
        .multiply(&b_thu_at_r)?
        .multiply_by_scalar(c_thu * determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)? * alpha);

    let b_thv_at_r = strain_displacement_matrix_thv_at_r(node_1_number, node_2_number, r, nodes)?;
    let b_thv_t_at_r = b_thv_at_r.transpose();
    let c_thv = young_modulus * i22_p;
    let k_thv_at_ip = b_thv_t_at_r
        .multiply(&b_thv_at_r)?
        .multiply_by_scalar(c_thv * determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)? * alpha);

    let b_thw_at_r = strain_displacement_matrix_thw_at_r(node_1_number, node_2_number, r, nodes)?;
    let b_thw_t_at_r = b_thw_at_r.transpose();
    let c_thw = young_modulus * i11_p;
    let k_thw_at_ip = b_thw_t_at_r
        .multiply(&b_thw_at_r)?
        .multiply_by_scalar(c_thw * determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)? * alpha);
    
    Ok(
        k_u_at_ip
            .add(&k_v_at_ip)?
            .add(&k_w_at_ip)?
            .add(&k_thu_at_ip)?
            .add(&k_thv_at_ip)?
            .add(&k_thw_at_ip)?
            .try_into_square_matrix()?
    )
}


fn compose_local_stiffness_matrix<V>(
    integration_points: &[(V, V)], 
    node_1_number: u32,
    node_2_number: u32,
    young_modulus: V,
    poisson_ratio: V,
    area: V,
    i11_p: V,
    i22_p: V,
    it: V,
    shear_factor: V,
    nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<SquareMatrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let mut local_stiffness_matrix = SquareMatrix::create(
        BEAM_NODES_NUMBER * BEAM_NODE_DOF, 
        &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
    );

    for (r, alpha) in integration_points
    {
        let local_stiffness_matrix_at_ip = local_stiffness_matrix_at_ip(
            node_1_number,
            node_2_number,
            young_modulus,
            poisson_ratio,
            area,
            i11_p,
            i22_p,
            it,
            shear_factor,
            *r,
            *alpha,
            nodes,
        )?;
        local_stiffness_matrix = local_stiffness_matrix.add(&local_stiffness_matrix_at_ip)?;
    }

    Ok(local_stiffness_matrix)
}


fn compose_rotation_matrix<V>(rotation_matrix_elements: &[V; 9]) -> SquareMatrix<V>
    where V: FloatTrait<>
{
    let [r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33] = 
        rotation_matrix_elements;

    SquareMatrix::create(
        BEAM_NODES_NUMBER * BEAM_NODE_DOF,
        &[
            [*r_11, *r_12, *r_13], [V::from(0f32); BEAM_NODE_DOF / 2], 
            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],

            [*r_21, *r_22, *r_23], [V::from(0f32); BEAM_NODE_DOF / 2], 
            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],

            [*r_31, *r_32, *r_33], [V::from(0f32); BEAM_NODE_DOF / 2],
            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],

            [V::from(0f32); BEAM_NODE_DOF / 2], [*r_11, *r_12, *r_13], 
            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],

            [V::from(0f32); BEAM_NODE_DOF / 2], [*r_21, *r_22, *r_23],
            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],

            [V::from(0f32); BEAM_NODE_DOF / 2], [*r_31, *r_32, *r_33],
            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],

            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
            [*r_11, *r_12, *r_13], [V::from(0f32); BEAM_NODE_DOF / 2],

            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
            [*r_21, *r_22, *r_23], [V::from(0f32); BEAM_NODE_DOF / 2],

            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
            [*r_31, *r_32, *r_33], [V::from(0f32); BEAM_NODE_DOF / 2],

            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
            [V::from(0f32); BEAM_NODE_DOF / 2], [*r_11, *r_12, *r_13],

            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
            [V::from(0f32); BEAM_NODE_DOF / 2], [*r_21, *r_22, *r_23],

            [V::from(0f32); BEAM_NODE_DOF / 2], [V::from(0f32); BEAM_NODE_DOF / 2],
            [V::from(0f32); BEAM_NODE_DOF / 2], [*r_31, *r_32, *r_33],
        ].concat(),
    )
}


pub struct Beam<V>
{
    node_1_number: u32,
    node_2_number: u32,
    young_modulus: V,
    poisson_ratio: V,
    area: V,
    i11_p: V,
    i22_p: V,
    it: V,
    shear_factor: V,
    rotation_matrix_elements: [V; 9],
    integration_points: [(V, V); 1],
}


impl<V> Beam<V>
    where V: FloatTrait<Output = V>
{
    pub fn create(
        number: u32,
        node_1_number: u32,
        node_2_number: u32,
        young_modulus: V,
        poisson_ratio: V,
        area: V,
        i11: V,
        i22: V,
        i12: V,
        it: V,
        shear_factor: V,
        local_axis_1_direction: [V; 3],
        nodes: &HashMap<u32, Node<V>>,
        rel_tol: V,
        abs_tol: V,
    )
        -> Result<Self, String>
    {
        check_beam_properties(
            number, 
            young_modulus, 
            poisson_ratio, 
            area, 
            i11, 
            i22, 
            it, 
            shear_factor, 
            node_1_number, 
            node_2_number, 
            nodes, 
            &local_axis_1_direction,
        )?;

        let (i11_p, i22_p, angle) = find_principal_moments_of_inertia(i11, i22, i12);

        let rotation_matrix_elements = find_rotation_matrix_elements(
            node_1_number, node_2_number, &local_axis_1_direction, angle, nodes, rel_tol, abs_tol,
        )?;
        let integration_points = [(V::from(0f32), V::from(2f32))];

        Ok(
            Beam 
            { 
                node_1_number, node_2_number, young_modulus, poisson_ratio, area, i11_p, i22_p, it, shear_factor, 
                rotation_matrix_elements, integration_points,
            }
        )
    }


    pub fn is_nodes_numbers_same(&self, node_1_number: u32, node_2_number: u32) -> bool
    {
        (node_1_number == self.node_1_number && node_2_number == self.node_2_number) ||
        (node_1_number == self.node_2_number && node_2_number == self.node_1_number)
    }


    pub fn extract_rotation_matrix(&self) -> SquareMatrix<V>
    {
        compose_rotation_matrix(&self.rotation_matrix_elements)
    }


    pub fn extract_local_stiffness_matrix(&self, nodes: &HashMap<u32, Node<V>>) -> Result<SquareMatrix<V>, String>
    {
        compose_local_stiffness_matrix(
            &self.integration_points,
            self.node_1_number,
            self.node_2_number,
            self.young_modulus,
            self.poisson_ratio,
            self.area,
            self.i11_p,
            self.i22_p,
            self.it,
            self.shear_factor,
            nodes,
        )
    }


    pub fn convert_uniformly_distributed_line_load_to_nodal_loads(
        &self,
        uniformly_distributed_line_load_value: V,
        nodes: &HashMap<u32, Node<V>>,
    ) 
        -> Result<Vector<V>, String>
    {
        let distributed_force_matrix = Matrix::create(
            1, 
            1,
            &[uniformly_distributed_line_load_value],
        );

        let mut nodal_forces = Vector::create(&[V::from(0f32); 2]);
        for (r, alpha) in self.integration_points.iter()
        {
            let determinant_of_jacobian_at_r = determinant_of_jacobian_at_r(
                self.node_1_number, self.node_2_number, *r, nodes,
            )?;

            let displacement_interpolation_matrix = Vector::create(&[h1_r(*r), h2_r(*r)]);

            let matrix = displacement_interpolation_matrix
                .multiply(&distributed_force_matrix)?
                .multiply_by_scalar(determinant_of_jacobian_at_r * *alpha);
            nodal_forces = nodal_forces.add(&matrix)?;
        }

        Ok(nodal_forces)
    }


    pub fn get_nodes_numbers(&self) -> [u32; 2]
    {
        [self.node_1_number, self.node_2_number]
    }


    pub fn extract_element_analysis_result(
        &self, nodes: &HashMap<u32, Node<V>>, displacements: &Vector<V>,
    )
        -> Result<Vec<(ElementForceComponent, V)>, String>
    {
        let node_1_index = nodes.get(&self.node_1_number)
            .ok_or(format!("Node: {} is absent!", self.node_1_number))?
            .get_index();
        let node_2_index = nodes.get(&self.node_2_number)
            .ok_or(format!("Node: {} is absent!", self.node_2_number))?
            .get_index();
        let mut global_displacements = Vector::create(
            &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
        );

        for i in 0..BEAM_NODE_DOF
        {
            *global_displacements.get_mut_element_value(&Position(i, 0))? = 
                *displacements.get_element_value(&Position(node_1_index * NODE_DOF + i, 0))?;
        }

        for i in 0..BEAM_NODE_DOF
        {
            *global_displacements.get_mut_element_value(&Position(i + BEAM_NODE_DOF, 0))? = 
                *displacements.get_element_value(&Position(node_2_index * NODE_DOF + i, 0))?;
        }

        let rotation_matrix = compose_rotation_matrix(&self.rotation_matrix_elements);
        let local_displacements = rotation_matrix.multiply(&global_displacements)?;

        let mut strain_displacement_matrix_u = Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
        );
        for (r, _alpha) in self.integration_points.iter()
        {
            let strain_displacement_matrix_u_at_r = strain_displacement_matrix_u_at_r(
                self.node_1_number, self.node_2_number, *r, nodes,
            )?;
            strain_displacement_matrix_u = strain_displacement_matrix_u.add(&strain_displacement_matrix_u_at_r)?;
        }
        let element_strains_u = strain_displacement_matrix_u.multiply(&local_displacements)?;
        let element_forces_u = element_strains_u.multiply_by_scalar(
            self.young_modulus * self.area / V::from(self.integration_points.len() as f32),
        );

        let shear_modulus = self.young_modulus / (V::from(2f32) * (V::from(1f32) + self.poisson_ratio));

        let mut strain_displacement_matrix_v = Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
        );
        for (r, _alpha) in self.integration_points.iter()
        {
            let strain_displacement_matrix_v_at_r = strain_displacement_matrix_v_at_r(
                self.node_1_number, self.node_2_number, *r, nodes,
            )?;
            strain_displacement_matrix_v = strain_displacement_matrix_v.add(&strain_displacement_matrix_v_at_r)?;
        }
        let element_strains_v = strain_displacement_matrix_v.multiply(&local_displacements)?;
        let element_forces_v = element_strains_v.multiply_by_scalar(
            shear_modulus * self.area  * self.shear_factor / V::from(self.integration_points.len() as f32),
        );
        let force_s_value = *element_forces_v.get_element_value(&Position(0, 0))?;

        let mut strain_displacement_matrix_w = Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
        );
        for (r, _alpha) in self.integration_points.iter()
        {
            let strain_displacement_matrix_w_at_r = strain_displacement_matrix_w_at_r(
                self.node_1_number, self.node_2_number, *r, nodes,
            )?;
            strain_displacement_matrix_w = strain_displacement_matrix_w.add(&strain_displacement_matrix_w_at_r)?;
        }
        let element_strains_w = strain_displacement_matrix_w.multiply(&local_displacements)?;
        let element_forces_w = element_strains_w.multiply_by_scalar(
            shear_modulus * self.area  * self.shear_factor / V::from(self.integration_points.len() as f32),
        );
        let force_t_value = *element_forces_w.get_element_value(&Position(0, 0))?;

        let mut strain_displacement_matrix_thu = Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
        );
        for (r, _alpha) in self.integration_points.iter()
        {
            let strain_displacement_matrix_thu_at_r = strain_displacement_matrix_thu_at_r(
                self.node_1_number, self.node_2_number, *r, nodes,
            )?;
            strain_displacement_matrix_thu = strain_displacement_matrix_thu.add(&strain_displacement_matrix_thu_at_r)?;
        }
        let element_strains_thu = strain_displacement_matrix_thu.multiply(&local_displacements)?;
        let element_forces_thu = element_strains_thu.multiply_by_scalar(
            shear_modulus * self.it / V::from(self.integration_points.len() as f32),
        );

        let beam_element_vector = find_2n_element_vector(self.node_1_number, self.node_2_number, nodes)?;
        let beam_element_length = beam_element_vector.norm()?;

        let mut strain_displacement_matrix_thv = Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
        );
        for (r, _alpha) in self.integration_points.iter()
        {
            let strain_displacement_matrix_thv_at_r = strain_displacement_matrix_thv_at_r(
                self.node_1_number, self.node_2_number, *r, nodes,
            )?;
            strain_displacement_matrix_thv = strain_displacement_matrix_thv.add(&strain_displacement_matrix_thv_at_r)?;
        }
        let element_strains_thv = strain_displacement_matrix_thv.multiply(&local_displacements)?;
        let element_forces_thv = element_strains_thv.multiply_by_scalar(
            self.young_modulus * self.i22_p / V::from(self.integration_points.len() as f32),
        );
        let moment_s_average = *element_forces_thv.get_element_value(&Position(0, 0))?;
        let moment_s_at_node_1 = moment_s_average + beam_element_length * force_t_value / V::from(2f32);
        let moment_s_at_node_2 = moment_s_average - beam_element_length * force_t_value / V::from(2f32);

        let mut strain_displacement_matrix_thw = Matrix::create(
            1,
            BEAM_NODES_NUMBER * BEAM_NODE_DOF,
            &[V::from(0f32); BEAM_NODES_NUMBER * BEAM_NODE_DOF],
        );
        for (r, _alpha) in self.integration_points.iter()
        {
            let strain_displacement_matrix_thw_at_r = strain_displacement_matrix_thw_at_r(
                self.node_1_number, self.node_2_number, *r, nodes,
            )?;
            strain_displacement_matrix_thw = strain_displacement_matrix_thw.add(&strain_displacement_matrix_thw_at_r)?;
        }
        let element_strains_thw = strain_displacement_matrix_thw.multiply(&local_displacements)?;
        let element_forces_thw = element_strains_thw.multiply_by_scalar(
            self.young_modulus * self.i11_p / V::from(self.integration_points.len() as f32),
        );
        let moment_t_average = *element_forces_thw.get_element_value(&Position(0, 0))?;
        let moment_t_at_node_1 = moment_t_average + beam_element_length * force_s_value / V::from(2f32);
        let moment_t_at_node_2 = moment_t_average - beam_element_length * force_s_value / V::from(2f32);

        let element_analysis_data = vec![
            (
                ElementForceComponent::ForceR,
                *element_forces_u.get_element_value(&Position(0, 0))?,
            ),
            (
                ElementForceComponent::ForceS,
                force_s_value,
            ),
            (
                ElementForceComponent::ForceT,
                force_t_value,
            ),
            (
                ElementForceComponent::MomentR,
                *element_forces_thu.get_element_value(&Position(0, 0))?,
            ),
            (
                ElementForceComponent::MomentS,
                moment_s_at_node_1,
            ),
            (
                ElementForceComponent::MomentS,
                moment_s_average,
            ),
            (
                ElementForceComponent::MomentS,
                moment_s_at_node_2,
            ),
            (
                ElementForceComponent::MomentT,
                moment_t_at_node_1,
            ),
            (
                ElementForceComponent::MomentT,
                moment_t_average,
            ),
            (
                ElementForceComponent::MomentT,
                moment_t_at_node_2,
            ),
        ];

        Ok(element_analysis_data)
    }
}
