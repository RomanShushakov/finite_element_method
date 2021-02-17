use crate::{FeNode, GlobalCoordinates};
use crate::{ElementsNumbers, ElementsValues};
use crate::TOLERANCE;
use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};
use crate::extended_matrix::extended_matrix::ExtendedMatrix;
use crate::extended_matrix::aux_traits_extended_matrix::{One};
use std::hash::Hash;
use std::fmt::Debug;


fn compare_with_tolerance(value: ElementsValues) -> ElementsValues
{
    if value.abs() < TOLERANCE { 0.0 } else { value }
}


struct TrussAuxFunctions<T, V>(T, V);


impl<T, V> TrussAuxFunctions<T, V>
    where T: Copy + From<ElementsNumbers> + Into<ElementsNumbers> + PartialOrd + Default +
             Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash +
             SubAssign + Debug + Mul<Output = T> + 'static,
          V: Copy + Into<ElementsValues> + From<ElementsValues> + Sub<Output = V> + Default +
             Mul<Output = V> + Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign +
             MulAssign + SubAssign + One + 'static
{
    fn length(node_1: &FeNode<T, V>, node_2: &FeNode<T, V>) -> V
    {
        V::from(((node_1.coordinates.x - node_2.coordinates.x).into().powi(2) +
        (node_1.coordinates.y - node_2.coordinates.y).into().powi(2) +
        (node_1.coordinates.z - node_2.coordinates.z).into().powi(2)).sqrt())
    }


    fn rotation_matrix(node_1: &FeNode<T, V>, node_2: &FeNode<T, V>) -> ExtendedMatrix<T, V>
    {
        let x = (node_2.coordinates.x - node_1.coordinates.x).into();
        let y = (node_2.coordinates.y - node_1.coordinates.y).into();
        let z = (node_2.coordinates.z - node_1.coordinates.z).into();
        let length = TrussAuxFunctions::<T, V>::length(node_1, node_2).into();
        let (u, v, w) = (length, 0.0, 0.0);
        let alpha = ((x * u + y * v + z * w) /
            (length * length)).acos();
        let (rotation_axis_coord_x, mut rotation_axis_coord_y,
            mut rotation_axis_coord_z) = (0.0 as ElementsValues, 0.0, 0.0);
        if x != 0.0 && y == 0.0 && z == 0.0
        {
            rotation_axis_coord_z = x;
        }
        else
        {
            rotation_axis_coord_y = z * length;
            rotation_axis_coord_z = - y * length;
        }
        let norm = 1.0 / (rotation_axis_coord_x.powi(2) +
            rotation_axis_coord_y.powi(2) + rotation_axis_coord_z.powi(2)).sqrt();
        let (x_n, y_n, z_n) = (rotation_axis_coord_x * norm,
            rotation_axis_coord_y * norm, rotation_axis_coord_z * norm);
        let (c, s) = (alpha.cos(), alpha.sin());
        let t = 1.0 - c;
        let q_11 = compare_with_tolerance(t * x_n * x_n + c);
        let q_12 = compare_with_tolerance(t * x_n * y_n - z_n * s);
        let q_13 = compare_with_tolerance(t * x_n * z_n + y_n * s);
        let q_21 = compare_with_tolerance(t * x_n * y_n + z_n * s);
        let q_22 = compare_with_tolerance(t * y_n * y_n + c);
        let q_23 = compare_with_tolerance(t * y_n * z_n - x_n * s);
        let q_31 = compare_with_tolerance(t * x_n * z_n - y_n * s);
        let q_32 = compare_with_tolerance(t * y_n * z_n + x_n * s);
        let q_33 = compare_with_tolerance(t * z_n * z_n + c);
        ExtendedMatrix::create(T::from(3 as ElementsNumbers),
           T::from(3 as ElementsNumbers),
           vec![V::from(q_11), V::from(q_12), V::from(q_13), V::from(q_21),
                V::from(q_22), V::from(q_23), V::from(q_31), V::from(q_32), V::from(q_33)])
    }


    fn power_func_x(a: V, x: V, n: i32) -> V
    {
        (0..n).fold(a, |acc, _| acc * x)
    }


    fn derivative_x(f: fn(V, V, i32) -> V,
                    a: V, x: V, n: i32) -> V
    {
        f(a * V::from(n as ElementsValues), x, n - 1)
    }


    fn dx_dr(x_1: V, x_2: V, r: V) -> V
    {
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5) * x_1,
            V::from(0.0), 0) -
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5) * x_1, r, 1) +
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5) * x_2,
            V::from(0.0), 0) +
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5) * x_2, r, 1)
    }


    fn jacobian(node_1: &FeNode<T, V>, node_2: &FeNode<T, V>, r: V) -> V
    {
        let length = TrussAuxFunctions::length(node_1, node_2);
        let x_1 = V::from(-1.0) * length / V::from(2.0);
        let x_2 = length / V::from(2.0);
        TrussAuxFunctions::<T, V>::dx_dr(x_1, x_2, r)
    }


    fn inverse_jacobian(node_1: &FeNode<T, V>, node_2: &FeNode<T, V>, r: V) -> V
    {
        V::from(1.0) / TrussAuxFunctions::jacobian(node_1, node_2, r)
    }


    fn determinant_of_jacobian(node_1: &FeNode<T, V>, node_2: &FeNode<T, V>, r: V) -> V
    {
        TrussAuxFunctions::jacobian(node_1, node_2, r)
    }


    fn dh1_dr(r: V) -> V
    {
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5), V::from(0.0), 0) -
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5), r, 1)
    }


    fn dh2_dr(r: V) -> V
    {
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5), V::from(0.0), 0) +
        TrussAuxFunctions::<T, V>::derivative_x(
            TrussAuxFunctions::<T, V>::power_func_x, V::from(0.5), r, 1)
    }
}


pub struct IntegrationPoint<V>
{
    sampling_point: V,
    weight: V,
}


pub struct State<T, V>
{
    pub rotation_matrix: ExtendedMatrix<T, V>,
    pub integration_points: Vec<IntegrationPoint<V>>
}


pub struct Truss2n2ip<'a, T, V>
{
    pub number: T,
    pub node_1: &'a FeNode<T, V>,
    pub node_2: &'a FeNode<T, V>,
    pub young_modulus: V,
    pub area: V,
    pub area_2: Option<V>,
    pub state: State<T, V>
}


impl<'a, T, V> Truss2n2ip<'a, T, V>
    where T: Copy + From<ElementsNumbers> + Into<ElementsNumbers> + PartialOrd + Default +
             Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash +
             SubAssign + Debug + Mul<Output = T> + 'static,
          V: Copy + Into<ElementsValues> + From<ElementsValues> + Sub<Output = V> + Default +
             Mul<Output = V> + Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign +
             MulAssign + SubAssign + One + 'static
{
    pub fn create(number: T, node_1: &'a FeNode<T, V>, node_2: &'a FeNode<T, V>,
        young_modulus: V, area: V, area_2: Option<V>) -> Self
    {
        let integration_point_1 = IntegrationPoint {
            sampling_point: V::from(- 1.0 / (3.0 as ElementsValues).sqrt()), weight: V::from(1.0) };
        let integration_point_2 = IntegrationPoint {
            sampling_point: V::from(1.0 / (3.0 as ElementsValues).sqrt()), weight: V::from(1.0) };
        let rotation_matrix = TrussAuxFunctions::rotation_matrix(node_1, node_2);
        let integration_points = vec![integration_point_1, integration_point_2];

        let state = State { rotation_matrix, integration_points };
        Truss2n2ip { number, node_1, node_2, young_modulus, area, area_2, state }
    }


    pub fn update(&mut self, node_1: &'a FeNode<T, V>, node_2: &'a FeNode<T, V>,
        young_modulus: V, area: V, area_2: Option<V>)
    {
        let rotation_matrix = TrussAuxFunctions::rotation_matrix(node_1, node_2);
        self.node_1 = node_1;
        self.node_2 = node_2;
        self.young_modulus = young_modulus;
        self.area = area;
        self.area_2 = area_2;
        self.state.rotation_matrix = rotation_matrix;
    }
}


fn area(area_1: ElementsValues, area_2: Option<ElementsValues>, r: ElementsValues)
    -> ElementsValues
{
    if let Some(area_2) = area_2
    {
        area_1.sqrt() + (area_1.sqrt() - area_2.sqrt()) * (r + 1.0) / 2.0
    }
    else
    {
        area_1
    }
}


fn power_func_x(a: ElementsValues, x: ElementsValues, n: i32) -> ElementsValues
{
    (0..n).fold(a, |acc, _| acc * x)
}


fn derivative_x(f: fn(ElementsValues, ElementsValues, i32) -> ElementsValues,
                a: ElementsValues, x: ElementsValues, n: i32) -> ElementsValues
{
    f(a * n as ElementsValues, x, n - 1)
}


fn dh1_dr(r: ElementsValues) -> ElementsValues
{
    derivative_x(power_func_x, 0.5, 0.0, 0) -
    derivative_x(power_func_x, 0.5, r, 1)
}


fn dh2_dr(r: ElementsValues) -> ElementsValues
{
    derivative_x(power_func_x, 0.5, 0.0, 0) +
    derivative_x(power_func_x, 0.5, r, 1)
}


fn dx_dr(x_1: ElementsValues, x_2: ElementsValues, r: ElementsValues) -> ElementsValues
{
    derivative_x(power_func_x, 0.5 * x_1, 0.0, 0) -
    derivative_x(power_func_x, 0.5 * x_1, r, 1) +
    derivative_x(power_func_x, 0.5 * x_2, 0.0, 0) +
    derivative_x(power_func_x, 0.5 * x_2, r, 1)
}


// fn jacobian(x_1: f32, x_2: f32, r: f32) -> f32
// {
//     derivative_x(power_func_x, 0.5 * x_2, r, 1) -
//     derivative_x(power_func_x, 0.5 * x_1, r, 1)
// }