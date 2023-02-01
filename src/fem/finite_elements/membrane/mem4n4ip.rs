use std::collections::HashMap;
use std::any::Any;

use extended_matrix::{Position, Matrix, FloatTrait, BasicOperationsTrait};

use crate::fem::finite_elements::finite_element::FiniteElementTrait;
use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::membrane::quad_full_mem_aux_functions::QuadFullMemAuxFunctions;

use crate::fem::global_analysis::fe_stiffness::{StiffnessGroup, StiffnessType};
use crate::fem::global_analysis::fe_dof_parameter_data::DOFParameterData;
use crate::fem::global_analysis::fe_global_analysis_result::Displacements;

use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::element_analysis::fe_element_analysis_result::
{
    ElementAnalysisData, ElementForces,
};

use crate::fem::finite_elements::membrane::consts::{MEMBRANE_NODE_DOF, MEM4N4IP_NODES_NUMBER};
use crate::fem::finite_elements::functions::extract_unique_elements_of_rotation_matrix;


struct IntegrationPoint<V>
{
    r: V,
    s: V,
    weight_r: V,
    weight_s: V,
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
        nodes_dof_parameters_global: Vec<DOFParameterData>
    ) 
        -> Self
    {
        State { rotation_matrix, integration_points, local_stiffness_matrix, nodes_dof_parameters_global }
    }
}


pub struct Mem4n4ip<V>
{
    node_1_number: u32,
    node_2_number: u32,
    node_3_number: u32,
    node_4_number: u32,
    young_modulus: V,
    poisson_ratio: V,
    thickness: V,
    state: State<V>,
}


impl<V> Mem4n4ip<V>
    where V: FloatTrait<Output = V>
{
    pub fn create(
        node_1_number: u32, 
        node_2_number: u32, 
        node_3_number: u32, 
        node_4_number: u32, 
        young_modulus: V, 
        poisson_ratio: V, 
        thickness: V,
        ref_nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    ) 
        -> Result<Self, String>
    {
        let integration_point_1 = IntegrationPoint {
            r: V::from(1f32 / 3f32).my_sqrt() * V::from(1f32), 
            s: V::from(1f32 / 3f32).my_sqrt() * V::from(1f32),
            weight_r: V::from(1f32), weight_s: V::from(1f32),
        };
        let integration_point_2 = IntegrationPoint {
            r: V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32), 
            s: V::from(1f32 / 3f32).my_sqrt() * V::from(1f32),
            weight_r: V::from(1f32), weight_s: V::from(1f32),
        };
        let integration_point_3 = IntegrationPoint {
            r: V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32), 
            s: V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32),
            weight_r: V::from(1f32), weight_s: V::from(1f32),
        };
        let integration_point_4 = IntegrationPoint {
            r: V::from(1f32 / 3f32).my_sqrt() * V::from(1f32), 
            s: V::from(1f32 / 3f32).my_sqrt() * V::from(-1f32),
            weight_r: V::from(1f32), weight_s: V::from(1f32),
        };

        let integration_points = vec![
            integration_point_1, integration_point_2, integration_point_3, integration_point_4
        ];

        let rotation_matrix = QuadFullMemAuxFunctions::rotation_matrix(
            node_2_number, node_3_number, node_4_number, ref_nodes, rel_tol, abs_tol,
        )?;

        let mut local_stiffness_matrix = Matrix::create(
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
            &[V::from(0f32); (MEM4N4IP_NODES_NUMBER * MEMBRANE_NODE_DOF).pow(2)],
        );

        for i in 0..integration_points.len()
        {
            let r = integration_points[i].r;
            let s = integration_points[i].s;
            let alpha = integration_points[i].weight_r * integration_points[i].weight_s;
            let matrix = QuadFullMemAuxFunctions::local_stiffness_matrix(
                node_1_number, 
                node_2_number, 
                node_3_number, 
                node_4_number, 
                young_modulus, 
                poisson_ratio, 
                thickness, 
                alpha, 
                r, 
                s, 
                &local_stiffness_matrix, 
                ref_nodes, 
                &rotation_matrix, 
                rel_tol,
            )?;
            local_stiffness_matrix = matrix;
        }

        let mut nodes_dof_parameters = QuadFullMemAuxFunctions::<V>::compose_node_dof_parameters(
            node_1_number,
        )?;
        let node_2_dof_parameters = QuadFullMemAuxFunctions::<V>::compose_node_dof_parameters(
            node_2_number,
        )?;
        let node_3_dof_parameters = QuadFullMemAuxFunctions::<V>::compose_node_dof_parameters(
            node_3_number,
        )?;
        let node_4_dof_parameters = QuadFullMemAuxFunctions::<V>::compose_node_dof_parameters(
            node_4_number,
        )?;

        nodes_dof_parameters.extend(node_2_dof_parameters);
        nodes_dof_parameters.extend(node_3_dof_parameters);
        nodes_dof_parameters.extend(node_4_dof_parameters);

        let state = State::create(rotation_matrix, integration_points,
            local_stiffness_matrix, nodes_dof_parameters);

        Ok(Mem4n4ip { 
            node_1_number,
            node_2_number, 
            node_3_number, 
            node_4_number, 
            young_modulus, 
            poisson_ratio, 
            thickness, 
            state 
        })
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

        let mut rows_number = self.state.nodes_dof_parameters_global.len();

        let element_global_displacements = Matrix::create(
            rows_number,
            1, 
            &element_global_displacements_values,
        );

        let element_local_displacements = self.state.rotation_matrix.multiply(&element_global_displacements)?;
        Ok(element_local_displacements)
    }
}


impl<V> FiniteElementTrait<V> for Mem4n4ip<V>
    where V: FloatTrait<Output = V>
{
    fn update(
        &mut self, 
        nodes_numbers: Vec<u32>, 
        properties: Vec<V>, 
        ref_nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    ) 
        -> Result<(), String>
    {
        let node_1_number = nodes_numbers[0];

        let node_2_number = nodes_numbers[1];

        let node_3_number = nodes_numbers[2];

        let node_4_number = nodes_numbers[3];

        let young_modulus = properties[0];

        let poisson_ratio = properties[1];

        let thickness = properties[2];

        let rotation_matrix = QuadFullMemAuxFunctions::rotation_matrix(
            node_2_number, node_3_number, node_4_number, ref_nodes, rel_tol, abs_tol)?;

        let mut local_stiffness_matrix = Matrix::create(
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
            &[V::from(0f32); (MEM4N4IP_NODES_NUMBER * MEMBRANE_NODE_DOF).pow(2)]
        );

        for i in 0..self.state.integration_points.len()
        {
            let r = self.state.integration_points[i].r;
            let s = self.state.integration_points[i].s;
            let alpha = self.state.integration_points[i].weight_r * self.state.integration_points[i].weight_s;
            let matrix = QuadFullMemAuxFunctions::local_stiffness_matrix(
                node_1_number, 
                node_2_number, 
                node_3_number, 
                node_4_number, 
                young_modulus, 
                poisson_ratio, 
                thickness, 
                alpha, 
                r, 
                s, 
                &local_stiffness_matrix, 
                ref_nodes, 
                &rotation_matrix, 
                rel_tol,
            )?;
            local_stiffness_matrix = matrix;
        }

        let mut nodes_dof_parameters = QuadFullMemAuxFunctions::<V>::compose_node_dof_parameters(
            node_1_number,
        )?;
        let node_2_dof_parameters = QuadFullMemAuxFunctions::<V>::compose_node_dof_parameters(
            node_2_number,
        )?;
        let node_3_dof_parameters = QuadFullMemAuxFunctions::<V>::compose_node_dof_parameters(
            node_3_number,
        )?;
        let node_4_dof_parameters = QuadFullMemAuxFunctions::<V>::compose_node_dof_parameters(
            node_4_number,
        )?;

        nodes_dof_parameters.extend(node_2_dof_parameters);
        nodes_dof_parameters.extend(node_3_dof_parameters);
        nodes_dof_parameters.extend(node_4_dof_parameters);

        self.node_1_number = node_1_number;
        self.node_2_number = node_2_number;
        self.node_3_number = node_3_number;
        self.node_4_number = node_4_number;
        self.young_modulus = young_modulus;
        self.poisson_ratio = poisson_ratio;
        self.thickness = thickness;
        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        self.state.nodes_dof_parameters_global = nodes_dof_parameters;

        Ok(())
    }


    fn extract_stiffness_matrix(&self) -> Result<Matrix<V>, String>
    {
        let mut interim_matrix = self.state.rotation_matrix.clone().transpose();
        if let Ok(matrix) =
        interim_matrix.multiply(&self.state.local_stiffness_matrix)
        {
            if let Ok(matrix) =
            matrix.multiply(&self.state.rotation_matrix)
            {
                return Ok(matrix);
            }
        }
        Err("Mem4n4ip: Stiffness matrix cannot be extracted!".to_string())
    }


    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup>
    {
        let (rows_number, columns_number) = (
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
        );

        let mut positions_kuu_1_1 = Vec::new();
        let mut positions_kuu_1_2 = Vec::new();
        let mut positions_kuu_1_3 = Vec::new();
        let mut positions_kuu_1_4 = Vec::new();

        let mut positions_kuu_2_1 = Vec::new();
        let mut positions_kuu_2_2 = Vec::new();
        let mut positions_kuu_2_3 = Vec::new();
        let mut positions_kuu_2_4 = Vec::new();

        let mut positions_kuu_3_1 = Vec::new();
        let mut positions_kuu_3_2 = Vec::new();
        let mut positions_kuu_3_3 = Vec::new();
        let mut positions_kuu_3_4 = Vec::new();

        let mut positions_kuu_4_1 = Vec::new();
        let mut positions_kuu_4_2 = Vec::new();
        let mut positions_kuu_4_3 = Vec::new();
        let mut positions_kuu_4_4 = Vec::new();

        let mut i = 0;
        while i < rows_number * columns_number
        {
            let position = Position(i / columns_number, i % columns_number);

            let row = i / columns_number;
            let column = i % columns_number;

            if row < QuadFullMemAuxFunctions::<V>::node_dof() && column < QuadFullMemAuxFunctions::<V>::node_dof()
            {
                positions_kuu_1_1.push(position);
            }
            else if row < QuadFullMemAuxFunctions::<V>::node_dof() &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() &&
                column < QuadFullMemAuxFunctions::<V>::node_dof() * 2
            {
                positions_kuu_1_2.push(position);
            }
            else if row < QuadFullMemAuxFunctions::<V>::node_dof() &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                column < QuadFullMemAuxFunctions::<V>::node_dof() * 3
            {
                positions_kuu_1_3.push(position);
            }
            else if row < QuadFullMemAuxFunctions::<V>::node_dof() &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() * 3
            {
                positions_kuu_1_4.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() &&
                row < QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                column < QuadFullMemAuxFunctions::<V>::node_dof()
            {
                positions_kuu_2_1.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() &&
                row < QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() &&
                column < QuadFullMemAuxFunctions::<V>::node_dof() * 2
            {
                positions_kuu_2_2.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() &&
                row < QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                column < QuadFullMemAuxFunctions::<V>::node_dof() * 3
            {
                positions_kuu_2_3.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() &&
                row < QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() * 3
            {
                positions_kuu_2_4.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                row < QuadFullMemAuxFunctions::<V>::node_dof() * 3 &&
                column < QuadFullMemAuxFunctions::<V>::node_dof()
            {
                positions_kuu_3_1.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                row < QuadFullMemAuxFunctions::<V>::node_dof() * 3 &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() &&
                column < QuadFullMemAuxFunctions::<V>::node_dof() * 2
            {
                positions_kuu_3_2.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                row < QuadFullMemAuxFunctions::<V>::node_dof() * 3 &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                column < QuadFullMemAuxFunctions::<V>::node_dof() * 3
            {
                positions_kuu_3_3.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                row < QuadFullMemAuxFunctions::<V>::node_dof() * 3 &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() * 3
            {
                positions_kuu_3_4.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() * 3 &&
                column < QuadFullMemAuxFunctions::<V>::node_dof()
            {
                positions_kuu_4_1.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() * 3 &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() &&
                column < QuadFullMemAuxFunctions::<V>::node_dof() * 2
            {
                positions_kuu_4_2.push(position);
            }
            else if row >= QuadFullMemAuxFunctions::<V>::node_dof() * 3 &&
                column >= QuadFullMemAuxFunctions::<V>::node_dof() * 2 &&
                column < QuadFullMemAuxFunctions::<V>::node_dof() * 3
            {
                positions_kuu_4_3.push(position);
            }
            else
            {
                positions_kuu_4_4.push(position);
            }
            i += 1;
        }

        vec![
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_1_number,
                number_2: self.node_1_number, 
                positions: positions_kuu_1_1, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_1_number,
                number_2: self.node_2_number, 
                positions: positions_kuu_1_2, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_1_number,
                number_2: self.node_3_number, 
                positions: positions_kuu_1_3, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_1_number,
                number_2: self.node_4_number, 
                positions: positions_kuu_1_4, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_2_number,
                number_2: self.node_1_number, 
                positions: positions_kuu_2_1, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_2_number,
                number_2: self.node_2_number, 
                positions: positions_kuu_2_2, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_2_number,
                number_2: self.node_3_number, 
                positions: positions_kuu_2_3, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_2_number,
                number_2: self.node_4_number, 
                positions: positions_kuu_2_4, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_3_number,
                number_2: self.node_1_number, 
                positions: positions_kuu_3_1,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_3_number,
                number_2: self.node_2_number, 
                positions: positions_kuu_3_2,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_3_number,
                number_2: self.node_3_number,
                positions: positions_kuu_3_3,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_3_number,
                number_2: self.node_4_number,
                positions: positions_kuu_3_4,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_4_number,
                number_2: self.node_1_number, 
                positions: positions_kuu_4_1,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_4_number,
                number_2: self.node_2_number, 
                positions: positions_kuu_4_2,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_4_number,
                number_2: self.node_3_number, 
                positions: positions_kuu_4_3,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_4_number,
                number_2: self.node_4_number, 
                positions: positions_kuu_4_4,
            },
        ]
    }


    fn is_node_belongs_to_element(&self, node_number: u32) -> bool
    {
        self.node_1_number == node_number || self.node_2_number == node_number ||
        self.node_3_number == node_number || self.node_4_number == node_number
    }


    fn refresh(
        &mut self, 
        ref_nodes: &HashMap<u32, FENode<V>>, 
        rel_tol: V, 
        abs_tol: V,
    ) 
        -> Result<(), String>
    {
        let rotation_matrix = QuadFullMemAuxFunctions::rotation_matrix(
            self.node_2_number, self.node_3_number, self.node_4_number, ref_nodes, rel_tol, abs_tol,
        )?;

        let mut local_stiffness_matrix = Matrix::create(
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
            QuadFullMemAuxFunctions::<V>::nodes_number() * QuadFullMemAuxFunctions::<V>::node_dof(),
            &[V::from(0f32); (MEM4N4IP_NODES_NUMBER * MEMBRANE_NODE_DOF).pow(2)],
        );

        for i in 0..self.state.integration_points.len()
        {
            let r = self.state.integration_points[i].r;
            let s = self.state.integration_points[i].s;
            let alpha = self.state.integration_points[i].weight_r * self.state.integration_points[i].weight_s;
            let matrix = QuadFullMemAuxFunctions::local_stiffness_matrix(
                self.node_1_number, 
                self.node_2_number, 
                self.node_3_number, 
                self.node_4_number, 
                self.young_modulus, 
                self.poisson_ratio, 
                self.thickness, 
                alpha, 
                r, 
                s, 
                &local_stiffness_matrix, 
                ref_nodes, 
                &rotation_matrix, 
                rel_tol,
            )?;
            local_stiffness_matrix = matrix;
        }

        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        Ok(())
    }


    fn is_nodes_numbers_same(&self, nodes_numbers: Vec<u32>) -> bool
    {
        nodes_numbers.iter().all(|node_number| self.is_node_belongs_to_element(*node_number))
    }


    fn extract_element_analysis_data(
        &self, 
        global_displacements: &Displacements<V>,
        ref_nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
    ) 
        -> Result<ElementAnalysisData<V>, String>
    {
        let element_local_displacements = self.extract_local_displacements(global_displacements)?;

        let c_matrix_multiplier = self.young_modulus / (V::from(1f32) - self.poisson_ratio.my_powi(2));
        let mut c_matrix = Matrix::create(
                3, 
                3, 
                &[
                    V::from(1f32), self.poisson_ratio, V::from(0f32),
                    self.poisson_ratio, V::from(1f32), V::from(0f32),
                    V::from(0f32), V::from(0f32), (V::from(1f32) - self.poisson_ratio) / V::from(2f32),
                ], 
            )
            .multiply_by_scalar(c_matrix_multiplier);

        // let mut strains_values = Vec::new();
        // let mut strains_components = Vec::new();

        // let mut stresses_values = Vec::new();
        // let mut stresses_components = Vec::new();

        let mut forces_values = Vec::new();
        let mut forces_components = Vec::new();

        let mut force_x = V::from(0f32);
        let mut force_y = V::from(0f32);
        let mut force_xy = V::from(0f32);

        let local_nodes_coordinates = vec![
            (V::from(1f32), V::from(1f32)), 
            (V::from(-1f32), V::from(1f32)),
            (V::from(-1f32), V::from(-1f32)), 
            (V::from(1f32), V::from(-1f32))
        ];

        for i in 0..local_nodes_coordinates.len()
        {
            let r = local_nodes_coordinates[i].0;
            let s = local_nodes_coordinates[i].1;

            let strain_displacement_matrix_at_node = QuadFullMemAuxFunctions::strain_displacement_matrix(
                self.node_1_number, 
                self.node_2_number, 
                self.node_3_number, 
                self.node_4_number, 
                r, 
                s, 
                ref_nodes, 
                &self.state.rotation_matrix, 
                rel_tol,
            )?;
            let strains_matrix_at_node = strain_displacement_matrix_at_node.multiply(&element_local_displacements)?;
            // let strains_at_node = QuadFullMemAuxFunctions::extract_column_matrix_values(&strains_matrix_at_node)?;
            let stresses_matrix_at_node = c_matrix.multiply(&strains_matrix_at_node)?;
            let stresses_at_node = QuadFullMemAuxFunctions::extract_column_matrix_values(
                &stresses_matrix_at_node,
            )?;
            for (j, k) in (0..3).into_iter().zip([0, 4, 1].into_iter())
            {
                // let stress_strain_component = 
                //     StressStrainComponent::iterator().nth(k).ok_or("Mem4n4ip: Unknown stress/strain component number!")?;
                // strains_values.push(strains_at_node[j]);
                // strains_components.push(*stress_strain_component);
                // stresses_values.push(stresses_at_node[j]);
                // stresses_components.push(*stress_strain_component);

                match k
                {
                    0 => force_x += stresses_at_node[j] * self.thickness,
                    4 => force_y += stresses_at_node[j] * self.thickness,
                    1 => force_xy += stresses_at_node[j] * self.thickness,
                    _ => (),
                }
            }
        }

        forces_values.push(force_x / V::from(4f32));
        forces_components.push(ForceComponent::MembraneForceX);
        forces_values.push(force_y / V::from(4f32));
        forces_components.push(ForceComponent::MembraneForceY);
        forces_values.push(force_xy / V::from(4f32));
        forces_components.push(ForceComponent::MembraneForceXY);
    
        // let element_strains = ElementStrains::create(strains_values, strains_components);
        // let element_stresses = ElementStresses::create(stresses_values, stresses_components);
        let element_forces = ElementForces::create(forces_values, forces_components);

        let element_analysis_data = ElementAnalysisData::create(
            None, 
            None,
            Some(element_forces),
            None,
        );
        Ok(element_analysis_data)
    }


    
    fn copy_nodes_numbers(&self) -> Vec<u32>
    {
        vec![self.node_1_number, self.node_2_number, self.node_3_number, self.node_4_number]
    }


    fn extract_unique_elements_of_rotation_matrix(&self) -> Result<Vec<V>, String>
    {
        extract_unique_elements_of_rotation_matrix(&self.state.rotation_matrix)
    }

    
    fn copy_properties(&self) -> Vec<V>
    {
        vec![self.young_modulus, self.poisson_ratio, self.thickness]
    }


    fn as_any(&self) -> &dyn Any
    {
        self
    }
}
