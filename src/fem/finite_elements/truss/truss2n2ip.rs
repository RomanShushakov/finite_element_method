use crate::{FeNode, GlobalCoordinates};
use crate::{ElementsNumbers, ElementsValues};
use crate::TOLERANCE;
use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};
use crate::extended_matrix::extended_matrix::ExtendedMatrix;
use crate::extended_matrix::aux_traits_extended_matrix::{One};
use std::hash::Hash;
use std::fmt::Debug;


pub struct Truss2n2ip<'a, T, V>
{
    pub number: T,
    pub node_1: &'a FeNode<T, V>,
    pub node_2: &'a FeNode<T, V>,
    pub young_modulus: V,
    pub area: V,
    pub area_2: Option<V>,
}


impl<'a, T, V> Truss2n2ip<'a, T, V>
    where T: Copy + From<ElementsNumbers> + Into<ElementsNumbers> + PartialOrd + Default +
             Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash +
             SubAssign + Debug + Mul<Output = T> + 'static,
          V: Copy + Into<ElementsValues> + From<ElementsValues> + Sub<Output = V> + Default +
             Mul<Output = V> + Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign +
             MulAssign + SubAssign + One + 'static
{
    pub fn length(&self) -> V
    {
        V::from(((self.node_1.coordinates.x - self.node_2.coordinates.x).into().powi(2) +
        (self.node_1.coordinates.y - self.node_2.coordinates.y).into().powi(2) +
        (self.node_1.coordinates.z - self.node_2.coordinates.z).into().powi(2)).sqrt())
    }


    pub fn rotation_matrix(&self) -> ExtendedMatrix<T, V>
    {
        let x = (self.node_2.coordinates.x - self.node_1.coordinates.x).into();
        let y = (self.node_2.coordinates.y - self.node_1.coordinates.y).into();
        let z = (self.node_2.coordinates.z - self.node_1.coordinates.z).into();
        let length = self.length().into();
        let (u, v, w) = (length, 0f64, 0f64);
        let alpha = ((x * u + y * v + z * w) /
            (length * length)).acos();
        let (rotation_axis_coord_x, mut rotation_axis_coord_y,
            mut rotation_axis_coord_z) = (0f64, 0f64, 0f64);
        if x != 0f64 && y == 0f64 && z == 0f64
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
        let q_11 = compare_with_tolerance(t * x_n * x_n + c, TOLERANCE);
        let q_12 = compare_with_tolerance(t * x_n * y_n - z_n * s, TOLERANCE);
        let q_13 = compare_with_tolerance(t * x_n * z_n + y_n * s, TOLERANCE);
        let q_21 = compare_with_tolerance(t * x_n * y_n + z_n * s, TOLERANCE);
        let q_22 = compare_with_tolerance(t * y_n * y_n + c, TOLERANCE);
        let q_23 = compare_with_tolerance(t * y_n * z_n - x_n * s, TOLERANCE);
        let q_31 = compare_with_tolerance(t * x_n * z_n - y_n * s, TOLERANCE);
        let q_32 = compare_with_tolerance(t * y_n * z_n + x_n * s, TOLERANCE);
        let q_33 = compare_with_tolerance(t * z_n * z_n + c, TOLERANCE);
        ExtendedMatrix::create(T::from(3u16), T::from(3u16),
           vec![V::from(q_11), V::from(q_12), V::from(q_13), V::from(q_21),
                V::from(q_22), V::from(q_23), V::from(q_31), V::from(q_32), V::from(q_33)])
    }
}


fn compare_with_tolerance(value: f64, tolerance: f64) -> f64
{
    if value.abs() < tolerance { 0.0 } else { value }
}