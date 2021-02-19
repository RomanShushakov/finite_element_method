use crate::fem::{GlobalCoordinates, FeNode, StiffnessGroup};
use crate::fem::StiffnessType;
use crate::fem::compare_with_tolerance;
use crate::extended_matrix::{ExtendedMatrix, One};
use crate::{ElementsNumbers, ElementsValues};

use std::rc::Rc;
use std::hash::Hash;
use std::fmt::Debug;
use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};


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
        ExtendedMatrix::create(T::from(6),
           T::from(6),
           vec![
               [V::from(q_11), V::from(q_12), V::from(q_13)], [V::from(0.0); 3],
               [V::from(q_21), V::from(q_22), V::from(q_23)], [V::from(0.0); 3],
               [V::from(q_31), V::from(q_32), V::from(q_33)], [V::from(0.0); 3],
               [V::from(0.0); 3], [V::from(q_11), V::from(q_12), V::from(q_13)],
               [V::from(0.0); 3], [V::from(q_21), V::from(q_22), V::from(q_23)],
               [V::from(0.0); 3], [V::from(q_31), V::from(q_32), V::from(q_33)],
           ].concat())
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


    fn strain_displacement_matrix(node_1: &FeNode<T, V>, node_2: &FeNode<T, V>, r: V)
        -> ExtendedMatrix<T, V>
    {
        let elements = vec![TrussAuxFunctions::<T, V>::dh1_dr(r), V::from(0.0),
            V::from(0.0), TrussAuxFunctions::<T, V>::dh2_dr(r), V::from(0.0), V::from(0.0)];
        let mut matrix = ExtendedMatrix::create(T::from(1),
            T::from(6), elements);
        let inverse_jacobian = TrussAuxFunctions::inverse_jacobian(node_1, node_2, r);
        matrix.multiply_by_number(inverse_jacobian);
        matrix
    }


    fn area(area_1: V, area_2: Option<V>, r: V) -> V
    {
        if let Some(area_2) = area_2
        {
            V::from((area_1.into().sqrt() + (area_2.into().sqrt() - area_1.into().sqrt()) *
                (r.into() + 1.0) / 2.0).powi(2))
        }
        else
        {
            area_1
        }
    }


    fn local_stiffness_matrix<'a>(node_1: &FeNode<T, V>, node_2: &FeNode<T, V>,
        young_modulus: V, area_1: V, area_2: Option<V>, alpha: V, r: V,
        local_stiffness_matrix: &ExtendedMatrix<T, V>)
        -> Result<ExtendedMatrix<T, V>, &'a str>
    {
        let current_area = TrussAuxFunctions::<T, V>::area(area_1, area_2, r);
        let mut lhs_matrix =
            TrussAuxFunctions::strain_displacement_matrix(node_1, node_2, r);
        lhs_matrix.transpose();
        lhs_matrix.multiply_by_number(young_modulus * current_area);
        let rhs_matrix =
            TrussAuxFunctions::strain_displacement_matrix(node_1, node_2, r);
        if let Ok(mut matrix) = lhs_matrix.multiply_by_matrix(&rhs_matrix)
        {
            matrix.multiply_by_number(
                TrussAuxFunctions::determinant_of_jacobian(node_1, node_2, r) * alpha);
            if let Ok(matrix) = local_stiffness_matrix.add(&matrix)
            {
                return Ok(matrix);
            }
        }
        Err("Truss2n2ip: Local stiffness matrix cannot be calculated!")
    }
}


struct IntegrationPoint<V>
{
    r: V,
    weight: V,
}


pub struct State<T, V>
{
    pub rotation_matrix: ExtendedMatrix<T, V>,
    integration_points: Vec<IntegrationPoint<V>>,
    pub local_stiffness_matrix: ExtendedMatrix<T, V>,
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
        young_modulus: V, area: V, area_2: Option<V>) -> Result<Self, &'a str>
    {
        let integration_point_1 = IntegrationPoint {
            r: V::from(- 1.0 / (3.0 as ElementsValues).sqrt()), weight: V::from(1.0) };
        let integration_point_2 = IntegrationPoint {
            r: V::from(1.0 / (3.0 as ElementsValues).sqrt()), weight: V::from(1.0) };
        let rotation_matrix = TrussAuxFunctions::rotation_matrix(node_1, node_2);
        let integration_points = vec![integration_point_1, integration_point_2];
        let mut local_stiffness_matrix =
            ExtendedMatrix::create(
                T::from(6), T::from(6), vec![V::from(0.0); 36]);
        for integration_point in &integration_points
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                node_1, node_2, young_modulus, area, area_2,
                integration_point.weight, integration_point.r,
                &local_stiffness_matrix)?;
            local_stiffness_matrix = matrix;
        }
        let state = State { rotation_matrix, integration_points, local_stiffness_matrix };
        Ok(Truss2n2ip { number, node_1, node_2, young_modulus, area, area_2, state })
    }


    pub fn update(&mut self, node_1: &'a FeNode<T, V>, node_2: &'a FeNode<T, V>,
        young_modulus: V, area: V, area_2: Option<V>) -> Result<(), &'a str>
    {
        let rotation_matrix = TrussAuxFunctions::rotation_matrix(node_1, node_2);
        let mut local_stiffness_matrix = ExtendedMatrix::create(
                T::from(6), T::from(6), vec![V::from(0.0); 36]);
        for integration_point in &self.state.integration_points
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                node_1, node_2, young_modulus, area, area_2,
                integration_point.weight, integration_point.r,
                &local_stiffness_matrix)?;
            local_stiffness_matrix = matrix;
        }
        self.node_1 = node_1;
        self.node_2 = node_2;
        self.young_modulus = young_modulus;
        self.area = area;
        self.area_2 = area_2;
        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        Ok(())
    }


    pub fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &'a str>
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
        Err("Truss2n2ip: Stiffness matrix cannot be extracted!")
    }


    pub fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>
    {
        let (rows_number, columns_number) = (T::from(6), T::from(6));
        let mut indexes_1_1 = Vec::new();
        let mut indexes_1_2 = Vec::new();
        let mut indexes_2_1 = Vec::new();
        let mut indexes_2_2 = Vec::new();
        for i in 0..(rows_number * columns_number).into()
        {
            let row = T::from(i) / columns_number;
            let column = T::from(i) % columns_number;
            if row < T::from(3) && column < T::from(3)
            {
                indexes_1_1.push(T::from(i));
            }
            else if row < T::from(3) && column >= T::from(3)
            {
                indexes_1_2.push(T::from(i));
            }
            else if row >= T::from(3) && column < T::from(3)
            {
                indexes_2_1.push(T::from(i));
            }
            else
            {
                indexes_2_2.push(T::from(i));
            }
        }
        vec![StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1.number,
                number_2: self.node_1.number, indexes: indexes_1_1, },
            StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1.number,
                number_2: self.node_2.number, indexes: indexes_1_2, },
            StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2.number,
                number_2: self.node_1.number, indexes: indexes_2_1 },
            StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2.number,
                number_2: self.node_2.number, indexes: indexes_2_2 }, ]
    }
}
