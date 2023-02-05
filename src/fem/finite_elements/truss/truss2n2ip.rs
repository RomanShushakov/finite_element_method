use std::collections::HashMap;
use std::any::Any;

use extended_matrix::{Position, Matrix, FloatTrait, BasicOperationsTrait};

use crate::fem::finite_elements::finite_element::FiniteElementTrait;
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

use crate::fem::finite_elements::truss::consts::{TRUSS_NODE_DOF, TRUSS2N2IP_NODES_NUMBER};
use crate::fem::finite_elements::functions::extract_unique_elements_of_rotation_matrix;

struct IntegrationPoint<V>
{
    r: V,
    weight: V,
}


struct State<V>
{
    rotation_matrix: Matrix<V>,
    integration_points: Vec<IntegrationPoint<V>>,
    local_stiffness_matrix: Matrix<V>,
    nodes_dof_parameters_global: Vec<DOFParameterData>,
}


impl<V> State<V>
{
    fn create(
        rotation_matrix: Matrix<V>, 
        integration_points: Vec<IntegrationPoint<V>>,
        local_stiffness_matrix: Matrix<V>,
        nodes_dof_parameters_global: Vec<DOFParameterData>,
    ) 
        -> Self
    {
        State { rotation_matrix, integration_points, local_stiffness_matrix,
        nodes_dof_parameters_global }
    }
}


pub struct Truss2n2ip<V>
{
    node_1_number: u32,
    node_2_number: u32,
    young_modulus: V,
    area: V,
    area_2: Option<V>,
    state: State<V>,
}


impl<V> Truss2n2ip<V>
    where V: FloatTrait<Output = V>
{
    pub fn create(
        node_1_number: u32, 
        node_2_number: u32, 
        young_modulus: V, 
        area: V, 
        area_2: Option<V>,
        nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    ) 
        -> Result<Self, String>
    {
        let integration_point_1 = IntegrationPoint {
            r: V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32), weight: V::from(1.0),
        };
        let integration_point_2 = IntegrationPoint {
            r: V::from(1f32 / 3f32).my_sqrt(), weight: V::from(1.0),
        };

        let rotation_matrix = TrussAuxFunctions::<V>::rotation_matrix(
            node_1_number, node_2_number, nodes, rel_tol, abs_tol,
        )?;

        let integration_points = vec![integration_point_1, integration_point_2];

        let mut local_stiffness_matrix = Matrix::create(
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            &[V::from(0f32); (TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF).pow(2)],
        );

        for integration_point in integration_points.iter()
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                node_1_number, 
                node_2_number, 
                young_modulus, 
                area, 
                area_2,
                integration_point.weight, 
                integration_point.r, 
                &local_stiffness_matrix,
                nodes,
            )?;
            local_stiffness_matrix = matrix;
        }

        let mut nodes_dof_parameters = TrussAuxFunctions::<V>::compose_node_dof_parameters(
            node_1_number,
        )?;

        let node_2_dof_parameters = TrussAuxFunctions::<V>::compose_node_dof_parameters(
            node_2_number,
        )?;

        nodes_dof_parameters.extend(node_2_dof_parameters);

        let state = State::create(rotation_matrix, integration_points,
            local_stiffness_matrix, nodes_dof_parameters);

        Ok(Truss2n2ip { node_1_number, node_2_number, young_modulus, area, area_2, state })
    }


    fn extract_local_displacements(&self, global_displacements: &Displacements<V>) -> Result<Matrix<V>, String>
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

        let rows_number = self.state.nodes_dof_parameters_global.len();

        let element_global_displacements = Matrix::create(
            rows_number, 1, &element_global_displacements_values,
        );

        let element_local_displacements = self.state.rotation_matrix.multiply(&element_global_displacements)?;
        Ok(element_local_displacements)
    }
}


impl<V> FiniteElementTrait<V> for Truss2n2ip<V>
    where V: FloatTrait<Output = V>
{
    fn update(
        &mut self, 
        nodes_numbers: Vec<u32>, 
        properties: Vec<V>, 
        nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    ) 
        -> Result<(), String>
    {
        let node_1_number = nodes_numbers[0];

        let node_2_number = nodes_numbers[1];

        let young_modulus = properties[0];

        let area = properties[1];

        let area_2 = if properties.len() == 3 { Some(properties[2]) } else { None };

        let rotation_matrix = TrussAuxFunctions::rotation_matrix(
            node_1_number, node_2_number, nodes, rel_tol, abs_tol,
        )?;

        let mut local_stiffness_matrix = Matrix::create(
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            &[V::from(0f32); (TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF).pow(2)],
        );

        for integration_point in self.state.integration_points.iter()
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                node_1_number, 
                node_2_number, 
                young_modulus, 
                area, 
                area_2,
                integration_point.weight, 
                integration_point.r, 
                &local_stiffness_matrix,
                nodes
            )?;
            local_stiffness_matrix = matrix;
        }

        let mut nodes_dof_parameters = TrussAuxFunctions::<V>::compose_node_dof_parameters(
            node_1_number,
        )?;

        let node_2_dof_parameters = TrussAuxFunctions::<V>::compose_node_dof_parameters(
            node_2_number,
        )?;

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


    fn extract_stiffness_matrix(&self) -> Result<Matrix<V>, String>
    {
        let interim_matrix = self.state.rotation_matrix.clone().transpose();
        if let Ok(matrix) = interim_matrix.multiply(&self.state.local_stiffness_matrix)
        {
            if let Ok(matrix) = matrix.multiply(&self.state.rotation_matrix)
            {
                return Ok(matrix);
            }
        }
        Err("Truss2n2ip: Stiffness matrix cannot be extracted!".to_string())
    }


    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup>
    {
        let (rows_number, columns_number) = (
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
        );

        let mut positions_1_1 = Vec::new();
        let mut positions_1_2 = Vec::new();
        let mut positions_2_1 = Vec::new();
        let mut positions_2_2 = Vec::new();

        let mut i = 0;
        while i < rows_number * columns_number
        {
            let position = Position(i / columns_number, i % columns_number);

            let row = i / columns_number;
            let column = i % columns_number;

            if row < TrussAuxFunctions::<V>::node_dof() && column < TrussAuxFunctions::<V>::node_dof()
            {
                positions_1_1.push(position);
            }
            else if row < TrussAuxFunctions::<V>::node_dof() && column >= TrussAuxFunctions::<V>::node_dof()
            {
                positions_1_2.push(position);
            }
            else if row >= TrussAuxFunctions::<V>::node_dof() && column < TrussAuxFunctions::<V>::node_dof()
            {
                positions_2_1.push(position);
            }
            else
            {
                positions_2_2.push(position);
            }
            i += 1;
        }

        vec![
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_1_number,
                number_2: self.node_1_number, 
                positions: positions_1_1, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_1_number,
                number_2: self.node_2_number, 
                positions: positions_1_2, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_2_number,
                number_2: self.node_1_number, 
                positions: positions_2_1 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_2_number,
                number_2: self.node_2_number, 
                positions: positions_2_2,
            },
        ]
    }


    fn is_node_belongs_to_element(&self, node_number: u32) -> bool
    {
        self.node_1_number == node_number || self.node_2_number == node_number
    }


    fn refresh(&mut self, nodes: &HashMap<u32, FENode<V>>, rel_tol: V, abs_tol: V) -> Result<(), String>
    {
        let rotation_matrix = TrussAuxFunctions::rotation_matrix(
            self.node_1_number, self.node_2_number, nodes, rel_tol, abs_tol,
        )?;

        let mut local_stiffness_matrix = Matrix::create(
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            TrussAuxFunctions::<V>::nodes_number() * TrussAuxFunctions::<V>::node_dof(),
            &[V::from(0f32); (TRUSS2N2IP_NODES_NUMBER * TRUSS_NODE_DOF).pow(2)],
        );

        for integration_point in self.state.integration_points.iter()
        {
            let matrix = TrussAuxFunctions::local_stiffness_matrix(
                self.node_1_number, 
                self.node_2_number, 
                self.young_modulus, 
                self.area,
                self.area_2, 
                integration_point.weight, 
                integration_point.r,
                &local_stiffness_matrix, 
                nodes,
            )?;
            local_stiffness_matrix = matrix;
        }

        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        Ok(())
    }


    fn is_nodes_numbers_same(&self, nodes_numbers: Vec<u32>) -> bool
    {
        (nodes_numbers[0] == self.node_1_number && nodes_numbers[1] == self.node_2_number) ||
        (nodes_numbers[0] == self.node_2_number && nodes_numbers[1] == self.node_1_number)
    }


    fn extract_element_analysis_data(
        &self, 
        global_displacements: &Displacements<V>,
        nodes: &HashMap<u32, FENode<V>>,
        _rel_tol: V,
    ) 
        -> Result<ElementAnalysisData<V>, String>
    {
        let element_local_displacements = self.extract_local_displacements(global_displacements)?;
        if self.area_2.is_some()
        {
            let mut forces_components = Vec::new();
            let mut forces_values = Vec::new();
            let mut axial_force = V::from(0f32);
            for ip in self.state.integration_points.iter()
            {
                let strain_displacement_matrix = TrussAuxFunctions::strain_displacement_matrix(
                    self.node_1_number, self.node_2_number, ip.r, nodes,
                )?;

                let strains_matrix = strain_displacement_matrix.multiply(&element_local_displacements)?;

                let area = TrussAuxFunctions::<V>::area(self.area, self.area_2, ip.r);

                let force_x = TrussAuxFunctions::extract_column_matrix_values(&strains_matrix)?[0] *
                    area * self.young_modulus;

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
                None, None, Some(element_forces), None,
            );
            Ok(element_analysis_data)
        }
        else
        {
            let r = self.state.integration_points[0].r;
            let mut forces_values = Vec::new();
            let mut forces_components = Vec::new();

            let strain_displacement_matrix = TrussAuxFunctions::strain_displacement_matrix(
                self.node_1_number, self.node_2_number, r, nodes,
            )?;
            let strains_matrix = strain_displacement_matrix.multiply(&element_local_displacements)?;
            let force_x = TrussAuxFunctions::extract_column_matrix_values(&strains_matrix)?[0] * 
                self.area * self.young_modulus;
            forces_components.push(ForceComponent::ForceX);
            forces_values.push(force_x);

            let element_forces = ElementForces::create(forces_values,
                forces_components);

            let element_analysis_data = ElementAnalysisData::create(
                None, None, Some(element_forces), None,
            );
            Ok(element_analysis_data)
        }
    }


    fn copy_nodes_numbers(&self) -> Vec<u32>
    {
        vec![self.node_1_number, self.node_2_number]
    }


    fn extract_unique_elements_of_rotation_matrix(&self) -> Result<Vec<V>, String>
    {
        extract_unique_elements_of_rotation_matrix(&self.state.rotation_matrix)
    }


    fn copy_properties(&self) -> Vec<V>
    {
        if let Some(area_2) = self.area_2
        {
            vec![self.young_modulus, self.area, area_2]
        }
        else
        {
            vec![self.young_modulus, self.area]
        }
    }


    fn as_any(&self) -> &dyn Any
    {
        self
    }
}
