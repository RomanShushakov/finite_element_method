use crate::fem::{FiniteElementTrait};
use crate::fem::{GlobalCoordinates, FeNode, StiffnessGroup, FEData, Displacement};
use crate::fem::{StiffnessType, GlobalForceDisplacementComponent};
use crate::fem::compare_with_tolerance;
use crate::extended_matrix::{ExtendedMatrix, MatrixElementPosition};
use crate::{ElementsNumbers, ElementsValues};

use std::rc::Rc;
use std::hash::Hash;
use std::fmt::Debug;
use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};
use std::cell::RefCell;


const TRUSS_NODE_DOF: ElementsNumbers = 3;
const TRUSS2N2IP_NODES_NUMBER: ElementsNumbers = 2;


struct TrussAuxFunctions<T, V>(T, V);


impl<T, V> TrussAuxFunctions<T, V>
    where T: Copy + From<ElementsNumbers> + Into<ElementsNumbers> + PartialOrd + Default +
             Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash +
             SubAssign + Debug + Mul<Output = T> + 'static,
          V: Copy + Into<ElementsValues> + From<ElementsValues> + Sub<Output = V> + Default +
             Mul<Output = V> + Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign +
             MulAssign + SubAssign + 'static
{
    fn length(node_1: Rc<RefCell<FeNode<T, V>>>, node_2: Rc<RefCell<FeNode<T, V>>>) -> V
    {
        V::from(((node_1.as_ref().borrow().coordinates.x - node_2.as_ref().borrow().coordinates.x)
            .into().powi(2) +
        (node_1.as_ref().borrow().coordinates.y - node_2.as_ref().borrow().coordinates.y)
            .into().powi(2) +
        (node_1.as_ref().borrow().coordinates.z - node_2.as_ref().borrow().coordinates.z)
            .into().powi(2)).sqrt())
    }


    fn rotation_matrix(node_1: Rc<RefCell<FeNode<T, V>>>, node_2: Rc<RefCell<FeNode<T, V>>>)
        -> ExtendedMatrix<T, V>
    {
        let x = (node_2.as_ref().borrow().coordinates.x -
            node_1.as_ref().borrow().coordinates.x).into();
        let y = (node_2.as_ref().borrow().coordinates.y -
            node_1.as_ref().borrow().coordinates.y).into();
        let z = (node_2.as_ref().borrow().coordinates.z -
            node_1.as_ref().borrow().coordinates.z).into();
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
        ExtendedMatrix::create(T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF),
           T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF),
           vec![
               [V::from(q_11), V::from(q_12), V::from(q_13)], [V::from(0.0); TRUSS_NODE_DOF as usize],
               [V::from(q_21), V::from(q_22), V::from(q_23)], [V::from(0.0); TRUSS_NODE_DOF as usize],
               [V::from(q_31), V::from(q_32), V::from(q_33)], [V::from(0.0); TRUSS_NODE_DOF as usize],
               [V::from(0.0); TRUSS_NODE_DOF as usize], [V::from(q_11), V::from(q_12), V::from(q_13)],
               [V::from(0.0); TRUSS_NODE_DOF as usize], [V::from(q_21), V::from(q_22), V::from(q_23)],
               [V::from(0.0); TRUSS_NODE_DOF as usize], [V::from(q_31), V::from(q_32), V::from(q_33)],
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


    fn jacobian(node_1: Rc<RefCell<FeNode<T, V>>>, node_2: Rc<RefCell<FeNode<T, V>>>, r: V) -> V
    {
        let length = TrussAuxFunctions::length(node_1, node_2);
        let x_1 = V::from(-1.0) * length / V::from(2.0);
        let x_2 = length / V::from(2.0);
        TrussAuxFunctions::<T, V>::dx_dr(x_1, x_2, r)
    }


    fn inverse_jacobian(node_1: Rc<RefCell<FeNode<T, V>>>, node_2: Rc<RefCell<FeNode<T, V>>>, r: V)
        -> V
    {
        V::from(1.0) / TrussAuxFunctions::jacobian(node_1, node_2, r)
    }


    fn determinant_of_jacobian(node_1: Rc<RefCell<FeNode<T, V>>>,
       node_2: Rc<RefCell<FeNode<T, V>>>, r: V) -> V
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


    fn strain_displacement_matrix(node_1: Rc<RefCell<FeNode<T, V>>>,
        node_2: Rc<RefCell<FeNode<T, V>>>, r: V) -> ExtendedMatrix<T, V>
    {
        let elements = vec![TrussAuxFunctions::<T, V>::dh1_dr(r), V::from(0.0),
            V::from(0.0), TrussAuxFunctions::<T, V>::dh2_dr(r), V::from(0.0), V::from(0.0)];
        let mut matrix = ExtendedMatrix::create(T::from(1),
            T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF), elements);
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


    fn local_stiffness_matrix(node_1: Rc<RefCell<FeNode<T, V>>>,
        node_2: Rc<RefCell<FeNode<T, V>>>, young_modulus: V, area_1: V, area_2: Option<V>,
        alpha: V, r: V, local_stiffness_matrix: &ExtendedMatrix<T, V>)
        -> Result<ExtendedMatrix<T, V>, String>
    {
        let current_area = TrussAuxFunctions::<T, V>::area(area_1, area_2, r);
        let mut lhs_matrix =
            TrussAuxFunctions::strain_displacement_matrix(
                Rc::clone(&node_1), Rc::clone(&node_2), r);
        lhs_matrix.transpose();
        lhs_matrix.multiply_by_number(young_modulus * current_area);
        let rhs_matrix =
            TrussAuxFunctions::strain_displacement_matrix(
                Rc::clone(&node_1), Rc::clone(&node_2), r);
        return match lhs_matrix.multiply_by_matrix(&rhs_matrix)
        {
            Ok(mut matrix) =>
                {
                    matrix.multiply_by_number(
                        TrussAuxFunctions::determinant_of_jacobian(
                            node_1, node_2, r) * alpha);
                    match local_stiffness_matrix.add_matrix(&matrix)
                    {
                        Ok(matrix) => Ok(matrix),
                        Err(e) =>
                            Err(format!("Truss2n2ip: Local stiffness matrix cannot be \
                                calculated! Reason: {}", e)),
                    }
                },
            Err(e) => Err(format!("Truss2n2ip: Local stiffness matrix cannot be \
                                calculated! Reason: {}", e)),
        }
    }


    fn compose_default_node_displacements<'a>(node_number: T)
        -> Result<Vec<Displacement<T, V>>, &'a str>
    {
        let mut nodes_displacements = Vec::new();
        for dof in 0..TRUSS_NODE_DOF
        {
            let displacement_component =
                GlobalForceDisplacementComponent::iterator().nth(dof as usize)
                    .ok_or("Truss2n2ip: Could not find displacement component!")?;
            let displacement = Displacement { node_number,
                component: *displacement_component,
                value: V::default() };
            nodes_displacements.push(displacement);
        }
        Ok(nodes_displacements)
    }
}


struct IntegrationPoint<V>
{
    r: V,
    weight: V,
}


struct State<T, V>
{
    rotation_matrix: ExtendedMatrix<T, V>,
    integration_points: Vec<IntegrationPoint<V>>,
    local_stiffness_matrix: ExtendedMatrix<T, V>,
    nodes_global_displacements: Vec<Displacement<T, V>>,
}


pub struct Truss2n2ip<T, V>
{
    number: T,
    node_1: Rc<RefCell<FeNode<T, V>>>,
    node_2: Rc<RefCell<FeNode<T, V>>>,
    young_modulus: V,
    area: V,
    area_2: Option<V>,
    state: State<T, V>
}


impl<T, V> Truss2n2ip<T, V>
    where T: Copy + From<ElementsNumbers> + Into<ElementsNumbers> + PartialOrd + Default +
             Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash +
             SubAssign + Debug + Mul<Output = T> + 'static,
          V: Copy + Into<ElementsValues> + From<ElementsValues> + Sub<Output = V> + Default +
             Mul<Output = V> + Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign +
             MulAssign + SubAssign + 'static
{
    pub fn create(number: T, node_1: Rc<RefCell<FeNode<T, V>>>,
        node_2: Rc<RefCell<FeNode<T, V>>>, young_modulus: V, area: V, area_2: Option<V>)
        -> Result<Self, String>
    {
        let integration_point_1 = IntegrationPoint {
            r: V::from(- 1.0 / (3.0 as ElementsValues).sqrt()), weight: V::from(1.0) };
        let integration_point_2 = IntegrationPoint {
            r: V::from(1.0 / (3.0 as ElementsValues).sqrt()), weight: V::from(1.0) };
        let rotation_matrix =
            TrussAuxFunctions::rotation_matrix(Rc::clone(&node_1),
                                               Rc::clone(&node_2));
        let integration_points = vec![integration_point_1, integration_point_2];
        let mut local_stiffness_matrix =
            ExtendedMatrix::create(
                T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF),
                T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF),
                vec![V::from(0.0); (TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF *
                    TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF) as usize ]);
        for integration_point in &integration_points
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                Rc::clone(&node_1), Rc::clone(&node_2), young_modulus,
                area, area_2, integration_point.weight, integration_point.r,
                &local_stiffness_matrix)?;
            local_stiffness_matrix = matrix;
        }
        let mut nodes_displacements =
            TrussAuxFunctions::compose_default_node_displacements(node_1.as_ref()
                .borrow().number)?;
        let node_2_displacements =
            TrussAuxFunctions::compose_default_node_displacements(node_2.as_ref()
                .borrow().number)?;
        nodes_displacements.extend(node_2_displacements);
        let state = State { rotation_matrix, integration_points, local_stiffness_matrix,
            nodes_global_displacements: nodes_displacements
        };
        Ok(Truss2n2ip { number, node_1, node_2, young_modulus, area, area_2, state })
    }


    pub fn update(&mut self, node_1: Rc<RefCell<FeNode<T, V>>>, node_2: Rc<RefCell<FeNode<T, V>>>,
        young_modulus: V, area: V, area_2: Option<V>) -> Result<(), String>
    {
        let rotation_matrix =
            TrussAuxFunctions::rotation_matrix(Rc::clone(&node_1),
                                               Rc::clone(&node_2));
        let mut local_stiffness_matrix = ExtendedMatrix::create(
                T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF),
                T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF),
                vec![V::from(0.0); (TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF *
                    TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF) as usize]);
        for integration_point in &self.state.integration_points
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                Rc::clone(&node_1), Rc::clone(&node_2), young_modulus,
                area, area_2, integration_point.weight, integration_point.r,
                &local_stiffness_matrix)?;
            local_stiffness_matrix = matrix;
        }
        let mut nodes_displacements =
            TrussAuxFunctions::compose_default_node_displacements(node_1.as_ref()
                .borrow().number)?;
        let node_2_displacements =
            TrussAuxFunctions::compose_default_node_displacements(node_2.as_ref()
                .borrow().number)?;
        nodes_displacements.extend(node_2_displacements);
        self.node_1 = node_1;
        self.node_2 = node_2;
        self.young_modulus = young_modulus;
        self.area = area;
        self.area_2 = area_2;
        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        self.state.nodes_global_displacements = nodes_displacements;
        Ok(())
    }


    pub fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
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
        let (rows_number, columns_number) =
            (T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF),
             T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF));
        let mut positions_1_1 = Vec::new();
        let mut positions_1_2 = Vec::new();
        let mut positions_2_1 = Vec::new();
        let mut positions_2_2 = Vec::new();
        for i in 0..(rows_number * columns_number).into()
        {
            let position = MatrixElementPosition { row: T::from(i) / columns_number,
                column: T::from(i) % columns_number };
            let row = T::from(i) / columns_number;
            let column = T::from(i) % columns_number;
            if row < T::from(TRUSS_NODE_DOF) && column < T::from(TRUSS_NODE_DOF)
            {
                positions_1_1.push(position);
            }
            else if row < T::from(TRUSS_NODE_DOF) && column >= T::from(TRUSS_NODE_DOF)
            {
                positions_1_2.push(position);
            }
            else if row >= T::from(TRUSS_NODE_DOF) && column < T::from(TRUSS_NODE_DOF)
            {
                positions_2_1.push(position);
            }
            else
            {
                positions_2_2.push(position);
            }
        }
        vec![StiffnessGroup { stiffness_type: StiffnessType::Kuu,
                number_1: self.node_1.as_ref().borrow().number,
                number_2: self.node_1.as_ref().borrow().number, positions: positions_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu,
                number_1: self.node_1.as_ref().borrow().number,
                number_2: self.node_2.as_ref().borrow().number, positions: positions_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu,
                number_1: self.node_2.as_ref().borrow().number,
                number_2: self.node_1.as_ref().borrow().number, positions: positions_2_1
             },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu,
                number_1: self.node_2.as_ref().borrow().number,
                number_2: self.node_2.as_ref().borrow().number, positions: positions_2_2
             }, ]
    }


    pub fn node_belong_element(&self, node_number: T) -> bool
    {
        self.node_1.as_ref().borrow().number == node_number ||
        self.node_2.as_ref().borrow().number == node_number
    }


    pub fn refresh(&mut self) -> Result<(), String>
    {
        let rotation_matrix =
            TrussAuxFunctions::rotation_matrix(Rc::clone(&self.node_1),
                                               Rc::clone(&self.node_2));
        let mut local_stiffness_matrix = ExtendedMatrix::create(
                T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF),
                T::from(TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF),
                vec![V::from(0.0); (TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF *
                    TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF) as usize]);
        for integration_point in &self.state.integration_points
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                Rc::clone(&self.node_1), Rc::clone(&self.node_2),
                self.young_modulus, self.area, self.area_2, integration_point.weight,
                integration_point.r, &local_stiffness_matrix)?;
            local_stiffness_matrix = matrix;
        }
        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        Ok(())
    }
}


impl<T, V> FiniteElementTrait<T, V> for Truss2n2ip<T, V>
    where T: Copy + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> +
             Mul<Output = T> + From<ElementsNumbers> + Into<ElementsNumbers> + Eq + Hash + Debug +
             SubAssign + PartialOrd + Default + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             Into<ElementsValues> + From<ElementsValues> + SubAssign + AddAssign + MulAssign +
             PartialEq + Debug + Default + 'static,
{
    fn update(&mut self, data: FEData<T, V>) -> Result<(), String>
    {
        if data.properties.len() == 3
        {
            self.update(
                Rc::clone(&data.nodes[0]), Rc::clone(&data.nodes[1]),
                data.properties[0], data.properties[1],
                Some(data.properties[2]))?;
        }
        else
        {
            self.update(
                Rc::clone(&data.nodes[0]), Rc::clone(&data.nodes[1]),
                data.properties[0], data.properties[1],
                None)?;
        }
        Ok(())
    }


    fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
    {
        self.extract_stiffness_matrix()
    }


    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>
    {
        self.extract_stiffness_groups()
    }


    fn node_belong_element(&self, node_number: T) -> bool
    {
        self.node_belong_element(node_number)
    }


    fn refresh(&mut self) -> Result<(), String>
    {
        self.refresh()?;
        Ok(())
    }


    fn number_same(&self, number: T) -> bool
    {
        self.number == number
    }


    fn nodes_numbers_same(&self, nodes_numbers: Vec<T>) -> bool
    {
        (nodes_numbers[0] == self.node_1.as_ref().borrow().number &&
        nodes_numbers[1] == self.node_2.as_ref().borrow().number) ||
        (nodes_numbers[0] == self.node_2.as_ref().borrow().number &&
        nodes_numbers[1] == self.node_1.as_ref().borrow().number)
    }
}