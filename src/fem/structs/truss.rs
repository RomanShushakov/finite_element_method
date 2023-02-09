use std::fmt::Debug;
use std::collections::HashMap;

use extended_matrix::
{
    FloatTrait, Vector3, VectorTrait, BasicOperationsTrait, Position, Matrix, TryIntoSquareMatrixTrait, SquareMatrix,
    Vector
};

use crate::fem::structs::{Node, NODE_DOF};
use crate::fem::methods_for_element_analysis::ElementForceComponent;


const TRUSS_NODES_NUMBER: usize = 2;
pub const TRUSS_NODE_DOF: usize = 3;


enum TrussDataError<V>
{
    YoungModulus(V),
    Area(V),
    Area2(V),
}


impl<V> TrussDataError<V>
    where V: Debug
{
    fn compose_error_message(&self) -> String
    {
        match self
        {
            TrussDataError::YoungModulus(value) => format!("Young's modulus {value:?} is less or equal to zero!"),
            TrussDataError::Area(value) => format!("Area {value:?} is less or equal to zero!"),
            TrussDataError::Area2(value) => format!("Area2 {value:?} is less or equal to zero!"),
        }
    }
}


fn check_truss_properties<V>(young_modulus: V,area: V, optional_area_2: Option<V>) -> Result<(), String>
    where V: Debug + PartialEq + PartialOrd + From<f32>
{
    if young_modulus <= V::from(0f32)
    {
        return Err(TrussDataError::<V>::YoungModulus(young_modulus).compose_error_message());
    }
    if area <= V::from(0f32)
    {
        return Err(TrussDataError::<V>::Area(area).compose_error_message());
    }
    if let Some(area_2) = optional_area_2
    {
        if area_2 <= V::from(0f32)
        {
            return Err(TrussDataError::<V>::Area2(area_2).compose_error_message());
        }
    }
    Ok(())
}


fn find_truss_element_vector<V>(
    node_1_number: u32, node_2_number: u32, nodes: &HashMap<u32, Node<V>>,
)
    -> Result<Vector3<V>, String>
    where V: FloatTrait<Output = V>
{
    let node_1 = nodes.get(&node_1_number).ok_or(format!("Node {node_1_number} does not exist!"))?;
    let node_2 = nodes.get(&node_2_number).ok_or(format!("Node {node_2_number} does not exist!"))?;

    let truss_element_vector_components: [V; 3] = node_2
        .get_coordinates()
        .iter()
        .zip(node_1.get_coordinates())
        .map(|(n, m)| *n - m)
        .collect::<Vec<V>>()
        .try_into()
        .map_err(|e| format!("{e:?} could not be converted to arr[3]"))?;

    Ok(Vector3::create(&truss_element_vector_components))
}


fn find_rotation_matrix_elements<V>(
    node_1_number: u32,
    node_2_number: u32,
    nodes: &HashMap<u32, Node<V>>,
    rel_tol: V,
    abs_tol: V,
)
    -> Result<[V; 9], String>
    where V: FloatTrait<Output = V>
{
    let truss_element_vector = find_truss_element_vector(node_1_number, node_2_number, nodes)?;
    let truss_element_length = truss_element_vector.norm()?;
    let direction_vector = Vector3::create(
        &[truss_element_length, V::from(0f32), V::from(0f32)],
    );

    let rotation_matrix = truss_element_vector
        .rotation_matrix_to_align_with_vector(&direction_vector, rel_tol, abs_tol)?;

    let mut rotation_matrix_elements = [V::from(0f32); 9];
    for i in 0..9
    {
        let element_value = *rotation_matrix.get_element_value(&Position(i / 3, i % 3))?;
        rotation_matrix_elements[i] = element_value;
    }

    Ok(rotation_matrix_elements)
}


fn power_func_x<V>(a: V, x: V, n: i32) -> V
    where V: FloatTrait<Output = V>
{
    (0..n).fold(a, |acc, _| acc * x)
}


fn derivative_x<V>(f: fn(V, V, i32) -> V, a: V, x: V, n: i32) -> V
    where V: FloatTrait<Output = V>
{
    let mut converted_n = V::from(0f32);
    (0..n).for_each(|_| converted_n += V::from(1f32));
    f(a * converted_n, x, n - 1)
}


fn dx_dr<V>(x_1: V, x_2: V, r: V) -> V
    where V: FloatTrait<Output = V>
{
    derivative_x(power_func_x, x_1 * V::from(0.5f32), V::from(0f32), 0) -
    derivative_x(power_func_x, x_1 * V::from(0.5f32), r, 1) +
    derivative_x(power_func_x, x_2 * V::from(0.5f32), V::from(0f32), 0) +
    derivative_x(power_func_x, x_2 * V::from(0.5f32), r, 1)
}


fn jacobian_at_r<V>(node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, Node<V>>) -> Result<V, String>
    where V: FloatTrait<Output = V>
{
    let truss_element_vector = find_truss_element_vector(node_1_number, node_2_number, nodes)?;
    let truss_element_length = truss_element_vector.norm()?;
    let x_1 = V::from(-1f32) * truss_element_length / V::from(2f32);
    let x_2 = truss_element_length / V::from(2f32);
    Ok(dx_dr(x_1, x_2, r))
}


fn inverse_jacobian_at_r<V>(
    node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<V, String>
    where V: FloatTrait<Output = V>
{
    Ok(V::from(1f32) / jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
}


fn determinant_of_jacobian_at_r<V>(
    node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<V, String>
    where V: FloatTrait<Output = V>
{
    jacobian_at_r(node_1_number, node_2_number, r, nodes)
}


fn dh1_dr<V>(r: V) -> V
    where V: FloatTrait<Output = V>
{
    derivative_x(power_func_x, V::from(0.5f32), V::from(0f32), 0) -
    derivative_x(power_func_x, V::from(0.5f32), r, 1)
}


fn dh2_dr<V>(r: V) -> V
    where V: FloatTrait<Output = V>
{
    derivative_x(power_func_x, V::from(0.5f32), V::from(0f32), 0) +
    derivative_x(power_func_x, V::from(0.5f32), r, 1)
}


fn strain_displacement_matrix_at_r<V>(
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
                TRUSS_NODES_NUMBER * TRUSS_NODE_DOF,
                &[dh1_dr(r), V::from(0f32), V::from(0f32), dh2_dr(r), V::from(0f32), V::from(0f32)],
            )
            .multiply_by_scalar(inverse_jacobian)
    )
}


fn area_at_r<V>(area: V, optional_area_2: Option<V>, r: V) -> V
    where V: FloatTrait<Output = V>
{
    if let Some(area_2) = optional_area_2
    {
        (area_2 - area) / V::from(2f32) * r + area - (area_2 - area) / V::from(2f32) * V::from(-1f32)
    }
    else
    {
        area
    }
}


fn local_stiffness_matrix_at_ip<V>(
    node_1_number: u32, 
    node_2_number: u32, 
    young_modulus: V, 
    area: V,
    optional_area_2: Option<V>,  
    r: V, 
    alpha: V,
    nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<SquareMatrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let b_at_r = strain_displacement_matrix_at_r(node_1_number, node_2_number, r, nodes)?;
    let b_t_at_r = b_at_r.transpose();
    let c_at_r = area_at_r(area, optional_area_2, r) * young_modulus;

    Ok(
        b_t_at_r
            .multiply_by_scalar(c_at_r)
            .multiply(&b_at_r)?
            .multiply_by_scalar(determinant_of_jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
            .multiply_by_scalar(alpha)
            .try_into_square_matrix()?
    )
}


fn compose_local_stiffness_matrix<V>(
    integration_points: &[(V, V)], 
    node_1_number: u32,
    node_2_number: u32,
    young_modulus: V,
    area: V,
    optional_area_2: Option<V>,
    nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<SquareMatrix<V>, String>
    where V: FloatTrait<Output = V>
{
    let mut local_stiffness_matrix = SquareMatrix::create(
        TRUSS_NODES_NUMBER * TRUSS_NODE_DOF, 
        &[V::from(0f32); TRUSS_NODES_NUMBER * TRUSS_NODE_DOF],
    );

    for (r, alpha) in integration_points
    {
        let local_stiffness_matrix_at_ip = local_stiffness_matrix_at_ip(
            node_1_number, node_2_number, young_modulus, area, optional_area_2, *r, *alpha, nodes,
        )?;
        local_stiffness_matrix = local_stiffness_matrix.add(&local_stiffness_matrix_at_ip)?;
    }

    Ok(local_stiffness_matrix)
}


fn compose_rotation_matrix<V>(rotation_matrix_elements: &[V; 9]) -> SquareMatrix<V>
    where V: FloatTrait
{
    let [q_11, q_12, q_13, q_21, q_22, q_23, q_31, q_32, q_33] = 
        rotation_matrix_elements;
    SquareMatrix::create(
        TRUSS_NODES_NUMBER * TRUSS_NODE_DOF,
        &vec![
            [*q_11, *q_12, *q_13], [V::from(0f32); TRUSS_NODE_DOF],
            [*q_21, *q_22, *q_23], [V::from(0f32); TRUSS_NODE_DOF],
            [*q_31, *q_32, *q_33], [V::from(0f32); TRUSS_NODE_DOF],
            [V::from(0f32); TRUSS_NODE_DOF], [*q_11, *q_12, *q_13],
            [V::from(0f32); TRUSS_NODE_DOF], [*q_21, *q_22, *q_23],
            [V::from(0f32); TRUSS_NODE_DOF], [*q_31, *q_32, *q_33],
        ].concat(),
    )
}


pub struct Truss<V>
{
    node_1_number: u32,
    node_2_number: u32,
    young_modulus: V,
    area: V,
    optional_area_2: Option<V>,
    rotation_matrix_elements: [V; 9],
    integration_points: [(V, V); 1],
}


impl<V> Truss<V>
    where V: FloatTrait<Output = V>
{
    pub fn create(
        node_1_number: u32,
        node_2_number: u32,
        young_modulus: V,
        area: V,
        optional_area_2: Option<V>,
        nodes: &HashMap<u32, Node<V>>,
        rel_tol: V,
        abs_tol: V,
    )
        -> Result<Self, String>
    {
        check_truss_properties(young_modulus, area, optional_area_2)?;

        let rotation_matrix_elements = find_rotation_matrix_elements(
            node_1_number, node_2_number, nodes, rel_tol, abs_tol,
        )?;
        let integration_points = [(V::from(0f32), V::from(2f32))];

        Ok(
            Truss
            {
                node_1_number, node_2_number, young_modulus, area, optional_area_2, rotation_matrix_elements,
                integration_points,
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
            self.area, 
            self.optional_area_2,
            nodes,
        )
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
            &[V::from(0f32); TRUSS_NODES_NUMBER * TRUSS_NODE_DOF],
        );

        for i in 0..TRUSS_NODE_DOF
        {
            *global_displacements.get_mut_element_value(&Position(i, 0))? = 
                *displacements.get_element_value(&Position(node_1_index * NODE_DOF + i, 0))?;
        }

        for i in 0..TRUSS_NODE_DOF
        {
            *global_displacements.get_mut_element_value(&Position(i + 3, 0))? = 
                *displacements.get_element_value(&Position(node_2_index * NODE_DOF + i, 0))?;
        }

        let rotation_matrix = compose_rotation_matrix(&self.rotation_matrix_elements);
        let local_displacements = rotation_matrix.multiply(&global_displacements)?;
        let mut strain_displacement_matrix = Matrix::create(
            1,
            TRUSS_NODES_NUMBER * TRUSS_NODE_DOF,
            &[V::from(0f32); TRUSS_NODES_NUMBER * TRUSS_NODE_DOF],
        );
        let mut area = V::from(0f32);
        for (r, _alpha) in self.integration_points.iter()
        {
            let strain_displacement_matrix_at_r = strain_displacement_matrix_at_r(
                self.node_1_number, self.node_2_number, *r, nodes,
            )?;
            strain_displacement_matrix = strain_displacement_matrix.add(&strain_displacement_matrix_at_r)?;
            area += area_at_r(self.area, self.optional_area_2, *r);
        }
        let element_strains = strain_displacement_matrix.multiply(&local_displacements)?;
        let element_forces = element_strains.multiply_by_scalar(self.young_modulus * area);

        let element_analysis_data = vec![
            (
                ElementForceComponent::ForceR,
                *element_forces.get_element_value(&Position(0, 0))?,
            ),
        ];

        Ok(element_analysis_data)
    }
}
