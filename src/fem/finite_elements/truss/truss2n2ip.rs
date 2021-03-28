use crate::fem::{FiniteElementTrait, Displacements, ForceComponent};
use crate::fem::
    {
        FENode, StiffnessGroup, FEData, DOFParameterData, ElementAnalysisData
    };
use crate::fem::{StiffnessType, GlobalDOFParameter, StressStrainComponent};
use crate::fem::compare_with_tolerance;
use crate::extended_matrix::{ExtendedMatrix, MatrixElementPosition, extract_element_value};
use crate::{ElementsNumbers, ElementsValues};

use std::rc::Rc;
use std::hash::Hash;
use std::fmt::Debug;
use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};
use std::cell::RefCell;


const TRUSS_NODE_DOF: ElementsNumbers = 3;
const TRUSS_STRESS_STRAIN_COMPONENTS_NUMBERS: [ElementsNumbers; 1] = [0];
const POINTS_NUMBER_FOR_TAPERED_TRUSS: ElementsNumbers = 5;
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
    fn length(node_1: Rc<RefCell<FENode<T, V>>>, node_2: Rc<RefCell<FENode<T, V>>>) -> V
    {
        V::from(((node_1.as_ref().borrow().coordinates.x - node_2.as_ref().borrow().coordinates.x)
            .into().powi(2) +
        (node_1.as_ref().borrow().coordinates.y - node_2.as_ref().borrow().coordinates.y)
            .into().powi(2) +
        (node_1.as_ref().borrow().coordinates.z - node_2.as_ref().borrow().coordinates.z)
            .into().powi(2)).sqrt())
    }


    fn rotation_matrix(node_1: Rc<RefCell<FENode<T, V>>>, node_2: Rc<RefCell<FENode<T, V>>>)
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


    fn jacobian(node_1: Rc<RefCell<FENode<T, V>>>, node_2: Rc<RefCell<FENode<T, V>>>, r: V) -> V
    {
        let length = TrussAuxFunctions::length(node_1, node_2);
        let x_1 = V::from(-1.0) * length / V::from(2.0);
        let x_2 = length / V::from(2.0);
        TrussAuxFunctions::<T, V>::dx_dr(x_1, x_2, r)
    }


    fn inverse_jacobian(node_1: Rc<RefCell<FENode<T, V>>>, node_2: Rc<RefCell<FENode<T, V>>>, r: V)
                        -> V
    {
        V::from(1.0) / TrussAuxFunctions::jacobian(node_1, node_2, r)
    }


    fn determinant_of_jacobian(node_1: Rc<RefCell<FENode<T, V>>>,
                               node_2: Rc<RefCell<FENode<T, V>>>, r: V) -> V
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


    fn strain_displacement_matrix(node_1: Rc<RefCell<FENode<T, V>>>,
                                  node_2: Rc<RefCell<FENode<T, V>>>, r: V) -> ExtendedMatrix<T, V>
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
            // V::from((area_1.into().sqrt() + (area_2.into().sqrt() - area_1.into().sqrt()) *
            //     (r.into() + 1.0) / 2.0).powi(2))
            (area_2 - area_1) / V::from(2.0) * r + area_1 -
                (area_2 - area_1) / V::from(2.0) * V::from(-1.0)
        }
        else
        {
            area_1
        }
    }


    fn local_stiffness_matrix(node_1: Rc<RefCell<FENode<T, V>>>,
                              node_2: Rc<RefCell<FENode<T, V>>>, young_modulus: V, area_1: V, area_2: Option<V>,
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


    fn compose_node_dof_parameters<'a>(node_number: T)
        -> Result<Vec<DOFParameterData<T>>, &'a str>
    {
        let mut node_dof_parameters = Vec::new();
        for dof in 0..TRUSS_NODE_DOF
        {
            let dof_parameter =
                GlobalDOFParameter::iterator().nth(dof as usize)
                    .ok_or("Truss2n2ip: Could not find dof parameter!")?;
            let dof_parameter = DOFParameterData { node_number,
                dof_parameter: *dof_parameter
            };
            node_dof_parameters.push(dof_parameter);
        }
        Ok(node_dof_parameters)
    }


    fn extract_column_matrix_values(column_matrix: &ExtendedMatrix<T, V>) -> Vec<V>
    {
        let mut values = Vec::new();
        let shape = column_matrix.get_shape();
        let all_values = column_matrix.extract_all_elements_values();
        for row in 0..shape.0.into()
        {
            for column in 0..shape.1.into()
            {
                let value = extract_element_value(T::from(row),
                    T::from(column), &all_values);
                values.push(value);
            }
        }
        values
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
    nodes_dof_parameters_global: Vec<DOFParameterData<T>>,
}


pub struct Truss2n2ip<T, V>
{
    number: T,
    node_1: Rc<RefCell<FENode<T, V>>>,
    node_2: Rc<RefCell<FENode<T, V>>>,
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
    pub fn create(number: T, node_1: Rc<RefCell<FENode<T, V>>>,
                  node_2: Rc<RefCell<FENode<T, V>>>, young_modulus: V, area: V, area_2: Option<V>)
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
        let mut nodes_dof_parameters =
            TrussAuxFunctions::<T, V>::compose_node_dof_parameters(
                node_1.as_ref().borrow().number)?;
        let node_2_dof_parameters =
            TrussAuxFunctions::<T, V>::compose_node_dof_parameters(
                node_2.as_ref().borrow().number)?;
        nodes_dof_parameters.extend(node_2_dof_parameters);
        let state = State { rotation_matrix, integration_points, local_stiffness_matrix,
            nodes_dof_parameters_global: nodes_dof_parameters
        };
        Ok(Truss2n2ip { number, node_1, node_2, young_modulus, area, area_2, state })
    }


    fn extract_local_displacements(&self, global_displacements: &Displacements<T, V>)
        -> Result<ExtendedMatrix<T, V>, String>
    {
        let mut element_global_displacements_values = Vec::new();
        for lhs_dof_parameter_data in &self.state.nodes_dof_parameters_global
        {
            if let Some(position) = global_displacements.dof_parameters_data
                .iter()
                .position(|rhs_dof_parameter_data|
                    rhs_dof_parameter_data == lhs_dof_parameter_data)
            {
                let displacement_value = global_displacements.displacements_values[position];
                element_global_displacements_values.push(displacement_value);
            }
            else
            {
                element_global_displacements_values.push(V::default());
            }
        }
        let element_global_displacements = ExtendedMatrix::create(
            T::from(self.state.nodes_dof_parameters_global.len() as ElementsNumbers),
            T::from(1),
            element_global_displacements_values);
        let element_local_displacements =
            self.state.rotation_matrix.multiply_by_matrix(&element_global_displacements)?;
        Ok(element_local_displacements)
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
        let node_1 = Rc::clone(&data.nodes[0]);
        let node_2 = Rc::clone(&data.nodes[1]);
        let young_modulus = data.properties[0];
        let area = data.properties[1];
        let area_2 =
            if data.properties.len() == 3 { Some(data.properties[2]) } else { None };
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
        let mut nodes_dof_parameters =
            TrussAuxFunctions::<T, V>::compose_node_dof_parameters(
                node_1.as_ref().borrow().number)?;
        let node_2_dof_parameters =
            TrussAuxFunctions::<T, V>::compose_node_dof_parameters(
                node_2.as_ref().borrow().number)?;
        nodes_dof_parameters.extend(node_2_dof_parameters);
        self.node_1 = node_1;
        self.node_2 = node_2;
        self.young_modulus = young_modulus;
        self.area = area;
        self.area_2 = area_2;
        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        self.state.nodes_dof_parameters_global = nodes_dof_parameters;
        Ok(())
    }


    fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
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


    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>
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


    fn node_belong_element(&self, node_number: T) -> bool
    {
        self.node_1.as_ref().borrow().number == node_number ||
        self.node_2.as_ref().borrow().number == node_number
    }


    fn refresh(&mut self) -> Result<(), String>
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


    fn extract_element_analysis_data(&self, global_displacements: &Displacements<T, V>)
        -> Result<ElementAnalysisData<T, V>, String>
    {
        let element_local_displacements =
            self.extract_local_displacements(global_displacements)?;
        if self.area_2.is_some()
        {
            let mut forces_components = Vec::new();
            let mut forces_values = Vec::new();
            let mut axial_force = V::default();
            for ip in &self.state.integration_points
            {
                let strain_displacement_matrix =
                    TrussAuxFunctions::strain_displacement_matrix(
                        Rc::clone(&self.node_1),
                        Rc::clone(&self.node_2), ip.r);
                let strains_matrix =
                    strain_displacement_matrix.multiply_by_matrix(&element_local_displacements)?;
                let stresses_matrix =
                    {
                        let mut matrix = strains_matrix.clone();
                        matrix.multiply_by_number(self.young_modulus);
                        matrix
                    };
                let stresses_values =
                    TrussAuxFunctions::extract_column_matrix_values(&stresses_matrix);
                for stress in &stresses_values
                {
                    let area = TrussAuxFunctions::<T, V>::area(
                        self.area, self.area_2, ip.r);
                    axial_force += *stress * area;
                }
            }
            let integral_axial_force =
                axial_force / V::from(self.state.integration_points.len() as ElementsValues);
            forces_components.push(ForceComponent::Axial);
            forces_values.push(integral_axial_force);
            let mut strains_components = Vec::new();
            let mut strains_values = Vec::new();
            let mut stresses_components = Vec::new();
            let mut stresses_values = Vec::new();
            for point_number in 0..POINTS_NUMBER_FOR_TAPERED_TRUSS
            {
                for component_number in &TRUSS_STRESS_STRAIN_COMPONENTS_NUMBERS
                {
                    let stress_strain_component =
                        StressStrainComponent::iterator().nth(*component_number as usize)
                            .ok_or(format!("Truss2n2ip: Stress strain component number {} \
                                could not be extracted", component_number))?;
                    strains_components.push(*stress_strain_component);
                    stresses_components.push(*stress_strain_component);
                }

                let r = V::from(-1.0) + V::from(2.0) /
                    V::from(POINTS_NUMBER_FOR_TAPERED_TRUSS as ElementsValues - 1.0) *
                    V::from(point_number as ElementsValues);
                let area = TrussAuxFunctions::<T, V>::area(self.area, self.area_2, r);
                let stress_value = integral_axial_force / area;
                stresses_values.push(stress_value);
                let strain_value = stress_value / self.young_modulus;
                strains_values.push(strain_value);
            }
            let element_analysis_data = ElementAnalysisData::create(
                self.number, strains_values, strains_components,
                stresses_values, stresses_components, forces_values, forces_components);
            Ok(element_analysis_data)
        }
        else
        {
            let r = self.state.integration_points[0].r;
            let mut strains_components = Vec::new();
            let mut stresses_components = Vec::new();
            let mut forces_components = Vec::new();
            let mut forces_values = Vec::new();
            let strain_displacement_matrix =
                TrussAuxFunctions::strain_displacement_matrix(
                Rc::clone(&self.node_1), Rc::clone(&self.node_2), r);
            let strains_matrix =
                strain_displacement_matrix.multiply_by_matrix(&element_local_displacements)?;
            let stresses_matrix =
                {
                    let mut matrix = strains_matrix.clone();
                    matrix.multiply_by_number(self.young_modulus);
                    matrix
                };
            for component_number in &TRUSS_STRESS_STRAIN_COMPONENTS_NUMBERS
            {
                let stress_strain_component =
                    StressStrainComponent::iterator().nth(*component_number as usize)
                        .ok_or(format!("Truss2n2ip: Stress strain component number {} could \
                            not be extracted", component_number))?;
                strains_components.push(*stress_strain_component);
                stresses_components.push(*stress_strain_component);
            }
            let strains_values =
                TrussAuxFunctions::extract_column_matrix_values(&strains_matrix);
            let stresses_values =
                TrussAuxFunctions::extract_column_matrix_values(&stresses_matrix);
            for stress in &stresses_values
            {
                let area = TrussAuxFunctions::<T, V>::area(self.area, self.area_2, r);
                let axial_force = *stress * area;
                forces_components.push(ForceComponent::Axial);
                forces_values.push(axial_force);
            }
            let element_analysis_data = ElementAnalysisData::create(
                self.number, strains_values, strains_components,
                stresses_values, stresses_components, forces_values, forces_components);
            Ok(element_analysis_data)
        }
    }


    fn extract_fe_number(&self) -> ElementsNumbers
    {
        self.number.into()
    }


    fn extract_nodes_numbers(&self) -> Vec<ElementsNumbers>
    {
        let mut numbers = Vec::new();
        let node_1_number = self.node_1.borrow().extract_number();
        let node_2_number = self.node_2.borrow().extract_number();
        numbers.push(node_1_number);
        numbers.push(node_2_number);
        numbers
    }


    fn extract_fe_properties(&self) -> Vec<ElementsValues>
    {
        let mut properties = Vec::new();
        properties.push(self.young_modulus.into());
        properties.push(self.area.into());
        if let Some(area) = self.area_2
        {
            properties.push(area.into());
        }
        properties
    }
}
