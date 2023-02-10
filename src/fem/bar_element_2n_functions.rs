use extended_matrix::{FloatTrait, Vector3, VectorTrait};

use std::collections::HashMap;

use crate::fem::structs::Node;
use crate::fem::math_functions::{power_func_x, derivative_x};


pub fn find_2n_element_vector<V>(
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
    let bar_2n_element_vector = find_2n_element_vector(node_1_number, node_2_number, nodes)?;
    let bar_2n_element_length = bar_2n_element_vector.norm()?;
    let x_1 = V::from(-1f32) * bar_2n_element_length / V::from(2f32);
    let x_2 = bar_2n_element_length / V::from(2f32);
    Ok(dx_dr(x_1, x_2, r))
}


pub fn inverse_jacobian_at_r<V>(
    node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<V, String>
    where V: FloatTrait<Output = V>
{
    Ok(V::from(1f32) / jacobian_at_r(node_1_number, node_2_number, r, nodes)?)
}


pub fn determinant_of_jacobian_at_r<V>(
    node_1_number: u32, node_2_number: u32, r: V, nodes: &HashMap<u32, Node<V>>,
) 
    -> Result<V, String>
    where V: FloatTrait<Output = V>
{
    jacobian_at_r(node_1_number, node_2_number, r, nodes)
}


pub fn h1_r<V>(r: V) -> V
    where V: FloatTrait<Output = V>
{
    V::from(0.5f32) * (V::from(1f32) - r)
}


pub fn h2_r<V>(r: V) -> V
    where V: FloatTrait<Output = V>
{
    V::from(0.5f32) * (V::from(1f32) + r)
}


pub fn dh1_dr<V>(r: V) -> V
    where V: FloatTrait<Output = V>
{
    derivative_x(power_func_x, V::from(0.5f32), V::from(0f32), 0) -
    derivative_x(power_func_x, V::from(0.5f32), r, 1)
}


pub fn dh2_dr<V>(r: V) -> V
    where V: FloatTrait<Output = V>
{
    derivative_x(power_func_x, V::from(0.5f32), V::from(0f32), 0) +
    derivative_x(power_func_x, V::from(0.5f32), r, 1)
}
