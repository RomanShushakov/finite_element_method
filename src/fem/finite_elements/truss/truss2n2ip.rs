use std::hash::Hash;
use std::fmt::Debug;
use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};

use extended_matrix::basic_matrix::basic_matrix::MatrixElementPosition;
use extended_matrix::extended_matrix::ExtendedMatrix;

use crate::fem::finite_elements::finite_element::{FiniteElementTrait, FEType};
use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::truss::truss_aux_functions::TrussAuxFunctions;

use crate::fem::global_analysis::fe_stiffness::{StiffnessGroup, StiffnessType};
use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData};
use crate::fem::global_analysis::fe_global_analysis_result::Displacements;

use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::element_analysis::fe_element_analysis_result::
{
    ElementAnalysisData, ElementForces
};

use crate::my_float::MyFloatTrait;

use crate::fem::finite_elements::truss::consts::{TRUSS_NODE_DOF, TRUSS2N2IP_NODES_NUMBER};
use std::collections::HashMap;


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


impl<T, V> State<T, V>
{
    fn create(rotation_matrix: ExtendedMatrix<T, V>, integration_points: Vec<IntegrationPoint<V>>,
        local_stiffness_matrix: ExtendedMatrix<T, V>,
        nodes_dof_parameters_global: Vec<DOFParameterData<T>>) -> Self
    {
        State { rotation_matrix, integration_points, local_stiffness_matrix,
        nodes_dof_parameters_global }
    }
}


pub struct Truss2n2ip<T, V>
{
    node_1_number: T,
    node_2_number: T,
    young_modulus: V,
    area: V,
    area_2: Option<V>,
    state: State<T, V>,
}


impl<T, V> Truss2n2ip<T, V>
    where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
             From<u8> + 'static,
          V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + From<f32> + Add<Output = V> +
             Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
             MyFloatTrait + PartialOrd + 'static
{
    pub fn create(node_1_number: T, node_2_number: T, young_modulus: V, area: V, area_2: Option<V>,
        tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<Self, String>
    {
        let integration_point_1 = IntegrationPoint {
            r: V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32), weight: V::from(1.0) };
        let integration_point_2 = IntegrationPoint {
            r: V::from(1f32 / 3f32).my_sqrt(), weight: V::from(1.0) };

        let rotation_matrix = TrussAuxFunctions::<T, V>::rotation_matrix(
            node_1_number, node_2_number, tolerance, nodes);

        let integration_points = vec![integration_point_1, integration_point_2];

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
            TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF).pow(2)], tolerance);

        for integration_point in &integration_points
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                node_1_number, node_2_number, young_modulus, area, area_2,
                integration_point.weight, integration_point.r, &local_stiffness_matrix,
                tolerance, nodes)?;
            local_stiffness_matrix = matrix;
        }

        let mut nodes_dof_parameters =
            TrussAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;

        let node_2_dof_parameters =
            TrussAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;

        nodes_dof_parameters.extend(node_2_dof_parameters);

        let state = State::create(rotation_matrix, integration_points,
            local_stiffness_matrix, nodes_dof_parameters);

        Ok(Truss2n2ip { node_1_number, node_2_number, young_modulus, area, area_2, state })
    }


    fn extract_local_displacements(&self, global_displacements: &Displacements<T, V>, tolerance: V)
        -> Result<ExtendedMatrix<T, V>, String>
    {
        let mut element_global_displacements_values = Vec::new();
        for lhs_dof_parameter_data in &self.state.nodes_dof_parameters_global
        {
            if let Some(position) = global_displacements.dof_parameters_data()
                .iter()
                .position(|rhs_dof_parameter_data|
                    rhs_dof_parameter_data == lhs_dof_parameter_data)
            {
                let displacement_value = global_displacements.displacements_values()[position];
                element_global_displacements_values.push(displacement_value);
            }
            else
            {
                element_global_displacements_values.push(V::from(0f32));
            }
        }

        let mut m = 0usize;
        let mut rows_number = T::from(0u8);
        while m < self.state.nodes_dof_parameters_global.len()
        {
            m += 1usize;
            rows_number += T::from(1u8);
        }

        let element_global_displacements = ExtendedMatrix::create(rows_number,
            T::from(1u8), element_global_displacements_values,
            tolerance);

        let element_local_displacements =
            self.state.rotation_matrix.multiply_by_matrix(&element_global_displacements)?;
        Ok(element_local_displacements)
    }
}


impl<T, V> FiniteElementTrait<T, V> for Truss2n2ip<T, V>
    where T: Copy + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> +
             Mul<Output = T> + Eq + Hash + Debug + SubAssign + PartialOrd + AddAssign +
             From<u8> + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             Into<f64> + SubAssign + AddAssign + MulAssign + PartialEq + Debug +
             MyFloatTrait + PartialOrd + From<f32> + 'static,
{
    fn update(&mut self, nodes_numbers: Vec<T>, properties: Vec<V>, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        let node_1_number = nodes_numbers[0];

        let node_2_number = nodes_numbers[1];

        let young_modulus = properties[0];

        let area = properties[1];

        let area_2 =
            if properties.len() == 3 { Some(properties[2]) } else { None };

        let rotation_matrix = TrussAuxFunctions::rotation_matrix(
            node_1_number, node_2_number, tolerance, nodes);

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
            TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF).pow(2)], tolerance);

        for integration_point in &self.state.integration_points
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                node_1_number, node_2_number, young_modulus, area, area_2,
                integration_point.weight, integration_point.r, &local_stiffness_matrix,
                tolerance, nodes)?;
            local_stiffness_matrix = matrix;
        }

        let mut nodes_dof_parameters =
            TrussAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;

        let node_2_dof_parameters =
            TrussAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;

        nodes_dof_parameters.extend(node_2_dof_parameters);

        self.node_1_number = node_1_number;
        self.node_2_number = node_2_number;
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
            (TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
             TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof());

        let mut positions_1_1 = Vec::new();
        let mut positions_1_2 = Vec::new();
        let mut positions_2_1 = Vec::new();
        let mut positions_2_2 = Vec::new();

        let mut i = T::from(0u8);
        while i < rows_number * columns_number
        {
            let position = MatrixElementPosition::create(
                i / columns_number, i % columns_number);

            let row = i / columns_number;
            let column = i % columns_number;

            if row < TrussAuxFunctions::<T, V>::node_dof() && column < TrussAuxFunctions::<T, V>::node_dof()
            {
                positions_1_1.push(position);
            }
            else if row < TrussAuxFunctions::<T, V>::node_dof() && column >= TrussAuxFunctions::<T, V>::node_dof()
            {
                positions_1_2.push(position);
            }
            else if row >= TrussAuxFunctions::<T, V>::node_dof() && column < TrussAuxFunctions::<T, V>::node_dof()
            {
                positions_2_1.push(position);
            }
            else
            {
                positions_2_2.push(position);
            }
            i += T::from(1u8);
        }

        vec![StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                number_2: self.node_1_number, positions: positions_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                number_2: self.node_2_number, positions: positions_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                number_2: self.node_1_number, positions: positions_2_1 },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                number_2: self.node_2_number, positions: positions_2_2 },
        ]
    }


    fn is_node_belongs_to_element(&self, node_number: T) -> bool
    {
        self.node_1_number == node_number || self.node_2_number == node_number
    }


    fn refresh(&mut self, tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        let rotation_matrix = TrussAuxFunctions::rotation_matrix(
            self.node_1_number, self.node_2_number, tolerance, nodes);

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
            TrussAuxFunctions::<T, V>::nodes_number() * TrussAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF).pow(2)], tolerance);

        for integration_point in self.state.integration_points.iter()
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                self.node_1_number, self.node_2_number, self.young_modulus, self.area,
                self.area_2, integration_point.weight, integration_point.r,
                &local_stiffness_matrix, tolerance, nodes)?;
            local_stiffness_matrix = matrix;
        }

        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        Ok(())
    }


    fn is_nodes_numbers_same(&self, nodes_numbers: Vec<T>) -> bool
    {
        (nodes_numbers[0] == self.node_1_number && nodes_numbers[1] == self.node_2_number) ||
        (nodes_numbers[0] == self.node_2_number && nodes_numbers[1] == self.node_1_number)
    }


    fn extract_element_analysis_data(&self, fe_type: FEType,
        global_displacements: &Displacements<T, V>, tolerance: V, nodes: &HashMap<T, FENode<V>>)
        -> Result<ElementAnalysisData<T, V>, String>
    {
        let element_local_displacements =
            self.extract_local_displacements(global_displacements, tolerance)?;
        if self.area_2.is_some()
        {
            let mut forces_components = Vec::new();
            let mut forces_values = Vec::new();
            let mut axial_force = V::from(0f32);
            for ip in self.state.integration_points.iter()
            {
                let strain_displacement_matrix =
                    TrussAuxFunctions::strain_displacement_matrix(self.node_1_number,
                        self.node_2_number, ip.r, tolerance, nodes);

                let strains_matrix =
                    strain_displacement_matrix.multiply_by_matrix(
                        &element_local_displacements)?;

                let area = TrussAuxFunctions::<T, V>::area(self.area, self.area_2, ip.r);

                let force_x = TrussAuxFunctions::extract_column_matrix_values(
                    &strains_matrix)[0] * area * self.young_modulus;

                axial_force += force_x;
            }

            let mut coeff = V::from(0f32);
            (0..self.state.integration_points.len()).for_each(|_| coeff += V::from(1f32));
            let integral_axial_force = axial_force / coeff;

            forces_components.push(ForceComponent::ForceX);
            forces_values.push(integral_axial_force);

            let element_forces = ElementForces::create(forces_values,
                forces_components);

            let element_analysis_data = ElementAnalysisData::create(
                fe_type, None, None, Some(element_forces), None);
            Ok(element_analysis_data)
        }
        else
        {
            let r = self.state.integration_points[0].r;
            let mut forces_values = Vec::new();
            let mut forces_components = Vec::new();

            let strain_displacement_matrix =
                TrussAuxFunctions::strain_displacement_matrix(
                    self.node_1_number, self.node_2_number, r, tolerance, nodes);
            let strains_matrix = strain_displacement_matrix.multiply_by_matrix(
                &element_local_displacements)?;
            let force_x = TrussAuxFunctions::extract_column_matrix_values(
                &strains_matrix)[0] * self.area * self.young_modulus;
            forces_components.push(ForceComponent::ForceX);
            forces_values.push(force_x);

            let element_forces = ElementForces::create(forces_values,
                forces_components);

            let element_analysis_data = ElementAnalysisData::create(
                fe_type, None, None, Some(element_forces), None);
            Ok(element_analysis_data)
        }
    }


    fn extract_nodes_numbers(&self) -> Vec<T>
    {
        let mut numbers = Vec::new();
        let node_1_number = self.node_1_number;
        let node_2_number = self.node_2_number;
        numbers.push(node_1_number);
        numbers.push(node_2_number);
        numbers
    }


    fn extract_fe_properties(&self) -> Vec<V>
    {
        let mut properties = Vec::new();
        properties.push(self.young_modulus);
        properties.push(self.area);
        if let Some(area) = self.area_2
        {
            properties.push(area);
        }
        properties
    }
}
