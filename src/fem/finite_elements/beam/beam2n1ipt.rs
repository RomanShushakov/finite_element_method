use std::collections::HashMap;
use std::any::Any;

use extended_matrix::{Position, Matrix, FloatTrait, BasicOperationsTrait};

use crate::fem::finite_elements::finite_element::FiniteElementTrait;
use crate::fem::finite_elements::fe_node::FENode;

use crate::fem::global_analysis::fe_stiffness::{StiffnessGroup, StiffnessType};
use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData};
use crate::fem::global_analysis::fe_global_analysis_result::Displacements;

use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::element_analysis::fe_element_analysis_result::{ElementAnalysisData, ElementForces, NodalForces};

use crate::fem::finite_elements::beam::beam_aux_functions::BeamAuxFunctions;

use crate::fem::finite_elements::beam::consts::{ BEAM_NODE_DOF, BEAM2N1IPT_NODES_NUMBER };
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
        State { rotation_matrix, integration_points, local_stiffness_matrix, nodes_dof_parameters_global }
    }
}


pub struct Beam2n1ipT<V>
{
    node_1_number: u32,
    node_2_number: u32,
    young_modulus: V,
    poisson_ratio: V,
    area: V,
    i11_init: V,
    i22_init: V,
    i12_init: V,
    shear_factor: V,
    it: V,
    local_axis_1_direction: [V; 3],
    state: State<V>,
}


impl<V> Beam2n1ipT<V>
    where V: FloatTrait<Output = V>
{
    pub fn create(
        node_1_number: u32, 
        node_2_number: u32, 
        young_modulus: V, 
        poisson_ratio: V, 
        area: V,
        i11_init: V, 
        i22_init: V, 
        i12_init: V, 
        it: V, 
        shear_factor: V,
        local_axis_1_direction: [V; 3], 
        nodes: &HashMap<u32, FENode<V>>,
        rel_tol: V,
        abs_tol: V,
    )
        -> Result<Self, String>
    {
        let (i11, i22, angle) = BeamAuxFunctions::<V>::find_principal_moments_of_inertia(
            i11_init, i22_init, i12_init,
        );

        let integration_point_1 = IntegrationPoint { r: V::from(0f32), weight: V::from(2f32) };

        let rotation_matrix = BeamAuxFunctions::<V>::rotation_matrix(
            node_1_number, node_2_number, &local_axis_1_direction, angle, nodes, rel_tol, abs_tol,
        )?;

        let integration_points = vec![integration_point_1];

        let mut local_stiffness_matrix = Matrix::create(
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &[V::from(0f32); (BEAM2N1IPT_NODES_NUMBER * BEAM_NODE_DOF).pow(2)],
        );

        for integration_point in &integration_points
        {
            let matrix = BeamAuxFunctions::local_stiffness_matrix(
                node_1_number, 
                node_2_number, 
                young_modulus, 
                poisson_ratio, 
                area, 
                shear_factor,
                i11, 
                i22, 
                it,
                integration_point.weight, 
                integration_point.r,
                nodes
            )?;

            local_stiffness_matrix = local_stiffness_matrix
                .add(&matrix)
                .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        let mut nodes_dof_parameters = BeamAuxFunctions::<V>::compose_node_dof_parameters(
            node_1_number,
        )?;

        let node_2_dof_parameters = BeamAuxFunctions::<V>::compose_node_dof_parameters(
            node_2_number,
        )?;

        nodes_dof_parameters.extend(node_2_dof_parameters);

        let state = State::create(rotation_matrix, integration_points,
            local_stiffness_matrix, nodes_dof_parameters);

        Ok(Beam2n1ipT { 
            node_1_number, 
            node_2_number, 
            young_modulus, 
            poisson_ratio,
            area, 
            i11_init, 
            i22_init, 
            i12_init, 
            shear_factor, 
            it, 
            local_axis_1_direction, 
            state,
        })
    }


    fn extract_local_displacements(
        &self, 
        global_displacements: &Displacements<V>, 
    )
        -> Result<Matrix<V>, String>
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

        let element_local_displacements =
            self.state.rotation_matrix.multiply(&element_global_displacements)?;
        Ok(element_local_displacements)
    }


    pub fn convert_uniformly_distributed_line_force_to_nodal_forces(
        &self,
        uniformly_distributed_line_force_value: V,
        ref_nodes: &HashMap<u32, FENode<V>>,
    ) 
        -> Result<Matrix<V>, String>
    {
        let distributed_force_matrix = Matrix::create(
            1, 
            1,
            &[uniformly_distributed_line_force_value],
        );

        let mut nodal_forces = Matrix::create(
            2,
            1,
            &[V::from(0f32); 2],
        );
        for integration_point in self.state.integration_points.iter()
        {
            let r = integration_point.r;
            let alpha = integration_point.weight;

            let determinant_of_jacobian = BeamAuxFunctions::<V>::determinant_of_jacobian(
                self.node_1_number, self.node_2_number, r, ref_nodes,
            );

            let mut displacement_interpolation_matrix = Matrix::create(
                    1,
                    2, 
                    &[BeamAuxFunctions::<V>::h1_r(r), BeamAuxFunctions::<V>::h2_r(r)],
                )
                .transpose();

            let mut matrix = displacement_interpolation_matrix
                .multiply(&distributed_force_matrix)?
                .multiply_by_scalar(determinant_of_jacobian * alpha);
            nodal_forces = nodal_forces.add(&matrix)?;
        }

        Ok(nodal_forces)
    }
}


impl<V> FiniteElementTrait<V> for Beam2n1ipT<V>
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

        let poisson_ratio = properties[1];

        let area = properties[2];

        let i11_init = properties[3];

        let i22_init = properties[4];

        let i12_init = properties[5];

        let (i11, i22, angle) = BeamAuxFunctions::<V>::find_principal_moments_of_inertia(
            i11_init, i22_init, i12_init,
        );

        let it = properties[6];

        let shear_factor = properties[7];

        let local_axis_1_direction = [properties[8], properties[9], properties[10]];

        let rotation_matrix = BeamAuxFunctions::rotation_matrix(
            node_1_number, node_2_number, &local_axis_1_direction, angle, nodes, rel_tol, abs_tol,
        )?;

        let mut local_stiffness_matrix = Matrix::create(
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &[V::from(0f32); (BEAM2N1IPT_NODES_NUMBER * BEAM_NODE_DOF).pow(2)],
        );

        for integration_point in &self.state.integration_points
        {
            let matrix = BeamAuxFunctions::local_stiffness_matrix(
                node_1_number, 
                node_2_number, 
                young_modulus, 
                poisson_ratio, 
                area, 
                shear_factor,
                i11, 
                i22,
                it,
                integration_point.weight, 
                integration_point.r,
                nodes,
            )?;
            local_stiffness_matrix = local_stiffness_matrix
                .add(&matrix)
                .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        let mut nodes_dof_parameters = BeamAuxFunctions::<V>::compose_node_dof_parameters(
            node_1_number,
        )?;

        let node_2_dof_parameters = BeamAuxFunctions::<V>::compose_node_dof_parameters(
            node_2_number,
        )?;

        nodes_dof_parameters.extend(node_2_dof_parameters);

        self.node_1_number = node_1_number;
        self.node_2_number = node_2_number;
        self.young_modulus = young_modulus;
        self.poisson_ratio = poisson_ratio;
        self.area = area;
        self.i11_init = i11_init;
        self.i22_init = i22_init;
        self.i12_init = i12_init;
        self.shear_factor = shear_factor;
        self.it = it;
        self.local_axis_1_direction = local_axis_1_direction;
        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        self.state.nodes_dof_parameters_global = nodes_dof_parameters;

        Ok(())
    }


    fn extract_stiffness_matrix(&self) -> Result<Matrix<V>, String>
    {
        let mut interim_matrix = self.state.rotation_matrix
            .clone()
            .transpose();

        if let Ok(matrix) = interim_matrix.multiply(&self.state.local_stiffness_matrix)
        {
            if let Ok(matrix) = matrix.multiply(&self.state.rotation_matrix)
            {
                return Ok(matrix);
            }
        }
        Err("Beam2n2ipT: Stiffness matrix cannot be extracted!".to_string())
    }


    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup>
    {
        let (rows_number, columns_number) = (
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
        );

        let mut positions_kuu_1_1 = Vec::new();
        let mut positions_kuth_1_1 = Vec::new();
        let mut positions_kthu_1_1 = Vec::new();
        let mut positions_kthth_1_1 = Vec::new();
        let mut positions_kuu_1_2 = Vec::new();
        let mut positions_kuth_1_2 = Vec::new();
        let mut positions_kthu_1_2 = Vec::new();
        let mut positions_kthth_1_2 = Vec::new();
        let mut positions_kuu_2_1 = Vec::new();
        let mut positions_kuth_2_1 = Vec::new();
        let mut positions_kthu_2_1 = Vec::new();
        let mut positions_kthth_2_1 = Vec::new();
        let mut positions_kuu_2_2 = Vec::new();
        let mut positions_kuth_2_2 = Vec::new();
        let mut positions_kthu_2_2 = Vec::new();
        let mut positions_kthth_2_2 = Vec::new();

        let mut i = 0;
        while i < rows_number * columns_number
        {
            let position = Position(i / columns_number, i % columns_number);

            let row = i / columns_number;
            let column = i % columns_number;

            if row < BeamAuxFunctions::<V>::node_dof() / 2 && column < BeamAuxFunctions::<V>::node_dof() / 2
            {
                positions_kuu_1_1.push(position);
            }
            else if row < BeamAuxFunctions::<V>::node_dof() / 2 &&
                column >= BeamAuxFunctions::<V>::node_dof() / 2 &&
                column < BeamAuxFunctions::<V>::node_dof()
            {
                positions_kuth_1_1.push(position);
            }
            else if row < BeamAuxFunctions::<V>::node_dof() / 2 &&
                column >= BeamAuxFunctions::<V>::node_dof() &&
                column < (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2)
            {
                positions_kuu_1_2.push(position);
            }
            else if row < BeamAuxFunctions::<V>::node_dof() / 2 &&
                column >= (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2)
            {
                positions_kuth_1_2.push(position);
            }
            else if row >= BeamAuxFunctions::<V>::node_dof() / 2 &&
                row < BeamAuxFunctions::<V>::node_dof() &&
                column < BeamAuxFunctions::<V>::node_dof() / 2
            {
                positions_kthu_1_1.push(position);
            }
            else if row >= BeamAuxFunctions::<V>::node_dof() / 2 &&
                row < BeamAuxFunctions::<V>::node_dof() &&
                column >= BeamAuxFunctions::<V>::node_dof() / 2 &&
                column < BeamAuxFunctions::<V>::node_dof()
            {
                positions_kthth_1_1.push(position);
            }
            else if row >= BeamAuxFunctions::<V>::node_dof() / 2 &&
                row < BeamAuxFunctions::<V>::node_dof() &&
                column >= BeamAuxFunctions::<V>::node_dof() &&
                column < (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2)
            {
                positions_kthu_1_2.push(position);
            }
            else if row >= BeamAuxFunctions::<V>::node_dof() / 2 &&
                row < BeamAuxFunctions::<V>::node_dof() &&
                column >= (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2)
            {
                positions_kthth_1_2.push(position);
            }
            else if row >= BeamAuxFunctions::<V>::node_dof() &&
                row < (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2) &&
                column < BeamAuxFunctions::<V>::node_dof() / 2
            {
                positions_kuu_2_1.push(position);
            }
            else if row >= BeamAuxFunctions::<V>::node_dof() &&
                row < (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2) &&
                column >= BeamAuxFunctions::<V>::node_dof() / 2 &&
                column < BeamAuxFunctions::<V>::node_dof()
            {
                positions_kuth_2_1.push(position);
            }
            else if row >= BeamAuxFunctions::<V>::node_dof() &&
                row < (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2) &&
                column >= BeamAuxFunctions::<V>::node_dof() &&
                column < (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2)
            {
                positions_kuu_2_2.push(position);
            }
            else if row >= BeamAuxFunctions::<V>::node_dof() &&
                row < (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2) &&
                column >= (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2)
            {
                positions_kuth_2_2.push(position);
            }
            else if row >= (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2) &&
                column < BeamAuxFunctions::<V>::node_dof() / 2
            {
                positions_kthu_2_1.push(position);
            }
            else if row >= (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2) &&
                column >= BeamAuxFunctions::<V>::node_dof() / 2 &&
                column < BeamAuxFunctions::<V>::node_dof()
            {
                positions_kthth_2_1.push(position);
            }
            else if row >= (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2) &&
                column >= BeamAuxFunctions::<V>::node_dof() &&
                column < (BeamAuxFunctions::<V>::node_dof() + BeamAuxFunctions::<V>::node_dof() / 2)
            {
                positions_kthu_2_2.push(position);
            }
            else
            {
                positions_kthth_2_2.push(position);
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
                stiffness_type: StiffnessType::Kuth, 
                number_1: self.node_1_number,
                number_2: self.node_1_number, 
                positions: positions_kuth_1_1, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_1_number,
                number_2: self.node_2_number,
                positions: positions_kuu_1_2,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuth, 
                number_1: self.node_1_number,
                number_2: self.node_2_number, 
                positions: positions_kuth_1_2, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kthu, 
                number_1: self.node_1_number,
                number_2: self.node_1_number, 
                positions: positions_kthu_1_1, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kthth, 
                number_1: self.node_1_number,
                number_2: self.node_1_number, 
                positions: positions_kthth_1_1, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kthu, 
                number_1: self.node_1_number,
                number_2: self.node_2_number, 
                positions: positions_kthu_1_2, 
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kthth, 
                number_1: self.node_1_number,
                number_2: self.node_2_number, 
                positions: positions_kthth_1_2,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu, 
                number_1: self.node_2_number,
                number_2: self.node_1_number, 
                positions: positions_kuu_2_1,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuth, 
                number_1: self.node_2_number,
                number_2: self.node_1_number, 
                positions: positions_kuth_2_1,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuu,
                number_1: self.node_2_number,
                number_2: self.node_2_number, 
                positions: positions_kuu_2_2,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kuth, 
                number_1: self.node_2_number,
                number_2: self.node_2_number,
                positions: positions_kuth_2_2,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kthu, 
                number_1: self.node_2_number,
                number_2: self.node_1_number, 
                positions: positions_kthu_2_1,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kthth, 
                number_1: self.node_2_number,
                number_2: self.node_1_number, 
                positions: positions_kthth_2_1,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kthu, 
                number_1: self.node_2_number,
                number_2: self.node_2_number, 
                positions: positions_kthu_2_2,
            },
            StiffnessGroup { 
                stiffness_type: StiffnessType::Kthth, 
                number_1: self.node_2_number,
                number_2: self.node_2_number, 
                positions: positions_kthth_2_2,
            },
        ]
    }


    fn is_node_belongs_to_element(&self, node_number: u32) -> bool
    {
        self.node_1_number == node_number || self.node_2_number == node_number
    }


    fn refresh(&mut self, nodes: &HashMap<u32, FENode<V>>, rel_tol: V, abs_tol: V) -> Result<(), String>
    {
        let (i11, i22, angle) = BeamAuxFunctions::<V>::find_principal_moments_of_inertia(
            self.i11_init, self.i22_init, self.i12_init,
        );

        let rotation_matrix = BeamAuxFunctions::rotation_matrix(
            self.node_1_number, self.node_2_number, &self.local_axis_1_direction, angle, nodes, rel_tol, abs_tol
        )?;

        let mut local_stiffness_matrix = Matrix::create(
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            BeamAuxFunctions::<V>::nodes_number() * BeamAuxFunctions::<V>::node_dof(),
            &[V::from(0f32); (BEAM2N1IPT_NODES_NUMBER * BEAM_NODE_DOF).pow(2)],
        );

        for integration_point in self.state.integration_points.iter()
        {
            let matrix = BeamAuxFunctions::local_stiffness_matrix(
                self.node_1_number, 
                self.node_2_number, 
                self.young_modulus, 
                self.poisson_ratio,
                self.area, 
                self.shear_factor, 
                i11,
                i22,
                self.it,
                integration_point.weight, 
                integration_point.r, 
                nodes
            )?;
            local_stiffness_matrix = local_stiffness_matrix
                .add(&matrix)
                .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
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
        rel_tol: V,
    ) 
        -> Result<ElementAnalysisData<V>, String>
    {
        let element_local_displacements = self.extract_local_displacements(global_displacements)?;

        let r = self.state.integration_points[0].r;
        let mut forces_values = Vec::new();
        let mut forces_components = Vec::new();

        let strain_displacement_matrix_u = BeamAuxFunctions::strain_displacement_matrix_u(
            self.node_1_number, self.node_2_number, r, nodes,
        )?;
        let strains_matrix_u = strain_displacement_matrix_u.multiply(&element_local_displacements)?;
        let force_x = BeamAuxFunctions::extract_column_matrix_values(&strains_matrix_u)?[0] * 
            self.area * self.young_modulus;
        forces_components.push(ForceComponent::ForceX);
        forces_values.push(force_x);

        let strain_displacement_matrix_v = BeamAuxFunctions::strain_displacement_matrix_v(
            self.node_1_number, self.node_2_number, r, nodes,
        )?;
        let strains_matrix_v = strain_displacement_matrix_v.multiply(&element_local_displacements)?;
        let shear_modulus = self.young_modulus / (V::from(2f32) * (V::from(1f32) + self.poisson_ratio));
        let force_y = BeamAuxFunctions::extract_column_matrix_values(&strains_matrix_v)?[0] * 
            shear_modulus * self.area * self.shear_factor;
        forces_components.push(ForceComponent::ForceY);
        forces_values.push(force_y);

        let strain_displacement_matrix_w = BeamAuxFunctions::strain_displacement_matrix_w(
            self.node_1_number, self.node_2_number, r, nodes,
        )?;
        let strains_matrix_w = strain_displacement_matrix_w.multiply(&element_local_displacements)?;
        let force_z = BeamAuxFunctions::extract_column_matrix_values(&strains_matrix_w)?[0] * 
            shear_modulus * self.area * self.shear_factor;
        forces_components.push(ForceComponent::ForceZ);
        forces_values.push(force_z);

        let strain_displacement_matrix_thu = BeamAuxFunctions::strain_displacement_matrix_thu(
            self.node_1_number, self.node_2_number, r, nodes,
        )?;
        let strains_matrix_thu = strain_displacement_matrix_thu.multiply(&element_local_displacements)?;
        let moment_x = BeamAuxFunctions::extract_column_matrix_values(&strains_matrix_thu)?[0] * 
            shear_modulus * self.it;
        forces_components.push(ForceComponent::MomentX);
        forces_values.push(moment_x);

        let length = BeamAuxFunctions::length(self.node_1_number, self.node_2_number, nodes);

        let mut forces_values_for_node_1 = Vec::new();
        let mut forces_components_for_node_1 = Vec::new();
        let mut forces_values_for_node_2 = Vec::new();
        let mut forces_components_for_node_2 = Vec::new();

        let (i11, i22, _angle) = BeamAuxFunctions::<V>::find_principal_moments_of_inertia(
            self.i11_init, self.i22_init, self.i12_init,
        );

        let strain_displacement_matrix_thv = BeamAuxFunctions::strain_displacement_matrix_thv(
            self.node_1_number, self.node_2_number, r, nodes,
        )?;
        let strains_matrix_thv = strain_displacement_matrix_thv.multiply(&element_local_displacements)?;
        let moment_y_average = BeamAuxFunctions::extract_column_matrix_values(&strains_matrix_thv)?[0] * 
            self.young_modulus * i22;
        forces_components.push(ForceComponent::MomentY);
        forces_values.push(moment_y_average);

        let moment_y_at_node_1 = moment_y_average + length * force_z / V::from(2f32);
        let moment_y_at_node_2 = moment_y_average - length * force_z / V::from(2f32);
        forces_components_for_node_1.push(ForceComponent::MomentY);
        forces_components_for_node_2.push(ForceComponent::MomentY);
        forces_values_for_node_1.push(moment_y_at_node_1);
        forces_values_for_node_2.push(moment_y_at_node_2);

        let strain_displacement_matrix_thw = BeamAuxFunctions::strain_displacement_matrix_thw(
            self.node_1_number, self.node_2_number, r, nodes,
        )?;
        let strains_matrix_thw = strain_displacement_matrix_thw.multiply(&element_local_displacements)?;
        let moment_z_average = BeamAuxFunctions::extract_column_matrix_values(&strains_matrix_thw)?[0] * 
            self.young_modulus * i11;
        forces_components.push(ForceComponent::MomentZ);
        forces_values.push(moment_z_average);

        let moment_z_at_node_1 = moment_z_average + length * force_y / V::from(2f32);
        let moment_z_at_node_2 = moment_z_average - length * force_y / V::from(2f32);
        forces_components_for_node_1.push(ForceComponent::MomentZ);
        forces_components_for_node_2.push(ForceComponent::MomentZ);
        forces_values_for_node_1.push(moment_z_at_node_1);
        forces_values_for_node_2.push(moment_z_at_node_2);

        let element_forces = ElementForces::create(forces_values, forces_components);

        let nodal_forces_for_node_1 = NodalForces::create(
            forces_values_for_node_1, forces_components_for_node_1,
        );
        let nodal_forces_for_node_2 = NodalForces::create(
            forces_values_for_node_2, forces_components_for_node_2,
        );
        let mut nodal_forces = HashMap::new();
        nodal_forces.insert(self.node_1_number, nodal_forces_for_node_1);
        nodal_forces.insert(self.node_2_number, nodal_forces_for_node_2);

        let element_analysis_data = ElementAnalysisData::create(
            None, 
            None, 
            Some(element_forces), 
            Some(nodal_forces),
        );
        Ok(element_analysis_data)

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
        vec![
            self.young_modulus, 
            self.poisson_ratio, 
            self.area, 
            self.i11_init, 
            self.i22_init,
            self.i12_init, 
            self.it, 
            self.shear_factor, 
            self.local_axis_1_direction[0],
            self.local_axis_1_direction[1], 
            self.local_axis_1_direction[2],
        ]
    }


    fn as_any(&self) -> &dyn Any
    {
        self
    }
}
