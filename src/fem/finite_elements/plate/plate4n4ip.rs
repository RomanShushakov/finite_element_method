use std::hash::Hash;
use std::fmt::Debug;
use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};
use std::collections::HashMap;

use extended_matrix::matrix_element_position::MatrixElementPosition;
use extended_matrix::extended_matrix::ExtendedMatrix;

use crate::fem::element_analysis::fe_stress_strain_components::StressStrainComponent;
use crate::fem::finite_elements::finite_element::{FiniteElementTrait, FEType};
use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::plate::quad_full_plate_aux_functions::QuadFullPlateAuxFunctions;

use crate::fem::global_analysis::fe_stiffness::{StiffnessGroup, StiffnessType};
use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData};
use crate::fem::global_analysis::fe_global_analysis_result::Displacements;

use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::element_analysis::fe_element_analysis_result::
{
    ElementAnalysisData, ElementForces, ElementStrains, ElementStresses
};

use extended_matrix_float::MyFloatTrait;

use crate::fem::finite_elements::plate::consts::{PLATE_NODE_DOF, PLATE4N4IP_NODES_NUMBER};
use crate::fem::finite_elements::functions::extract_unique_elements_of_rotation_matrix;


struct IntegrationPoint<V>
{
    r: V,
    s: V,
    weight_r: V,
    weight_s: V,
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


pub struct Plate4n4ip<T, V>
{
    node_1_number: T,
    node_2_number: T,
    node_3_number: T,
    node_4_number: T,
    young_modulus: V,
    poisson_ratio: V,
    thickness: V,
    shear_factor: V,
    state: State<T, V>,
}


impl<T, V> Plate4n4ip<T, V>
    where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
             From<u8> + Ord + 'static,
          V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + From<f32> + Add<Output = V> +
             Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
             MyFloatTrait + PartialOrd + 'static
{
    pub fn create(node_1_number: T, node_2_number: T, node_3_number: T, node_4_number: T, 
        young_modulus: V, poisson_ratio: V, thickness: V, shear_factor: V, tolerance: V, 
        ref_nodes: &HashMap<T, FENode<V>>) -> Result<Self, String>
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

        let rotation_matrix = QuadFullPlateAuxFunctions::rotation_matrix(
            node_2_number, node_3_number, node_4_number, tolerance, ref_nodes)?;

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(),
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (PLATE4N4IP_NODES_NUMBER * PLATE_NODE_DOF).pow(2)], tolerance)?;

        for i in 0..integration_points.len()
        {
            let r = integration_points[i].r;
            let s = integration_points[i].s;
            let alpha = integration_points[i].weight_r * integration_points[i].weight_s;

            let matrix = QuadFullPlateAuxFunctions::local_stiffness_matrix(
                node_1_number, node_2_number, node_3_number, node_4_number, young_modulus, poisson_ratio, 
                thickness, shear_factor, alpha, r, s, ref_nodes, &rotation_matrix, tolerance)?;

            local_stiffness_matrix = local_stiffness_matrix.add_matrix(&matrix)
                .map_err(|e| format!("Plate4n4ip: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        let mut nodes_dof_parameters =
            QuadFullPlateAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;
        let node_2_dof_parameters =
            QuadFullPlateAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;
        let node_3_dof_parameters =
            QuadFullPlateAuxFunctions::<T, V>::compose_node_dof_parameters(node_3_number)?;
        let node_4_dof_parameters =
            QuadFullPlateAuxFunctions::<T, V>::compose_node_dof_parameters(node_4_number)?;

        nodes_dof_parameters.extend(node_2_dof_parameters);
        nodes_dof_parameters.extend(node_3_dof_parameters);
        nodes_dof_parameters.extend(node_4_dof_parameters);

        let state = State::create(rotation_matrix, integration_points,
            local_stiffness_matrix, nodes_dof_parameters);

        Ok(Plate4n4ip { node_1_number, node_2_number, node_3_number, node_4_number,
            young_modulus, poisson_ratio, thickness, shear_factor, state })
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
            tolerance)?;

        let element_local_displacements =
            self.state.rotation_matrix.multiply_by_matrix(&element_global_displacements)?;
        Ok(element_local_displacements)
    }
}


impl<T, V> FiniteElementTrait<T, V> for Plate4n4ip<T, V>
    where T: Copy + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> +
             Mul<Output = T> + Eq + Hash + Debug + SubAssign + PartialOrd + AddAssign +
             From<u8> + Ord + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             Into<f64> + SubAssign + AddAssign + MulAssign + PartialEq + Debug +
             MyFloatTrait + PartialOrd + From<f32> + 'static,
{
    fn update(&mut self, nodes_numbers: Vec<T>, properties: Vec<V>, tolerance: V,
        ref_nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        let node_1_number = nodes_numbers[0];

        let node_2_number = nodes_numbers[1];

        let node_3_number = nodes_numbers[2];

        let node_4_number = nodes_numbers[3];

        let young_modulus = properties[0];

        let poisson_ratio = properties[1];

        let thickness = properties[2];

        let shear_factor = properties[3];

        let rotation_matrix = QuadFullPlateAuxFunctions::rotation_matrix(
            node_2_number, node_3_number, node_4_number, tolerance, ref_nodes)?;

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(),
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (PLATE4N4IP_NODES_NUMBER * PLATE_NODE_DOF).pow(2)], tolerance)?;

        for i in 0..self.state.integration_points.len()
        {
            let r = self.state.integration_points[i].r;
            let s = self.state.integration_points[i].s;
            let alpha = self.state.integration_points[i].weight_r * self.state.integration_points[i].weight_s;

            let matrix = QuadFullPlateAuxFunctions::local_stiffness_matrix(
                node_1_number, node_2_number, node_3_number, node_4_number, young_modulus, poisson_ratio, 
                thickness, shear_factor, alpha, r, s, ref_nodes, &rotation_matrix, tolerance)?;

            local_stiffness_matrix = local_stiffness_matrix.add_matrix(&matrix)
                .map_err(|e| format!("Plate4n4ip: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        let mut nodes_dof_parameters =
            QuadFullPlateAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;
        let node_2_dof_parameters =
            QuadFullPlateAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;
        let node_3_dof_parameters =
            QuadFullPlateAuxFunctions::<T, V>::compose_node_dof_parameters(node_3_number)?;
        let node_4_dof_parameters =
            QuadFullPlateAuxFunctions::<T, V>::compose_node_dof_parameters(node_4_number)?;

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
        self.shear_factor = shear_factor;
        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        self.state.nodes_dof_parameters_global = nodes_dof_parameters;

        Ok(())
    }


    fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, String>
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
        Err("Plate4n4ip: Stiffness matrix cannot be extracted!".to_string())
    }


    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>
    {
        let (rows_number, columns_number) =
            (QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(),
             QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof());

            let mut positions_kuu_1_1 = Vec::new();
            let mut positions_kuth_1_1 = Vec::new();
            let mut positions_kuu_1_2 = Vec::new();
            let mut positions_kuth_1_2 = Vec::new();
            let mut positions_kuu_1_3 = Vec::new();
            let mut positions_kuth_1_3 = Vec::new();
            let mut positions_kuu_1_4 = Vec::new();
            let mut positions_kuth_1_4 = Vec::new();

            let mut positions_kthu_1_1 = Vec::new();
            let mut positions_kthth_1_1 = Vec::new();
            let mut positions_kthu_1_2 = Vec::new();
            let mut positions_kthth_1_2 = Vec::new();
            let mut positions_kthu_1_3 = Vec::new();
            let mut positions_kthth_1_3 = Vec::new();
            let mut positions_kthu_1_4 = Vec::new();
            let mut positions_kthth_1_4 = Vec::new();

            let mut positions_kuu_2_1 = Vec::new();
            let mut positions_kuth_2_1 = Vec::new();
            let mut positions_kuu_2_2 = Vec::new();
            let mut positions_kuth_2_2 = Vec::new();
            let mut positions_kuu_2_3 = Vec::new();
            let mut positions_kuth_2_3 = Vec::new();
            let mut positions_kuu_2_4 = Vec::new();
            let mut positions_kuth_2_4 = Vec::new();

            let mut positions_kthu_2_1 = Vec::new();
            let mut positions_kthth_2_1 = Vec::new();
            let mut positions_kthu_2_2 = Vec::new();
            let mut positions_kthth_2_2 = Vec::new();
            let mut positions_kthu_2_3 = Vec::new();
            let mut positions_kthth_2_3 = Vec::new();
            let mut positions_kthu_2_4 = Vec::new();
            let mut positions_kthth_2_4 = Vec::new();

            let mut positions_kuu_3_1 = Vec::new();
            let mut positions_kuth_3_1 = Vec::new();
            let mut positions_kuu_3_2 = Vec::new();
            let mut positions_kuth_3_2 = Vec::new();
            let mut positions_kuu_3_3 = Vec::new();
            let mut positions_kuth_3_3 = Vec::new();
            let mut positions_kuu_3_4 = Vec::new();
            let mut positions_kuth_3_4 = Vec::new();

            let mut positions_kthu_3_1 = Vec::new();
            let mut positions_kthth_3_1 = Vec::new();
            let mut positions_kthu_3_2 = Vec::new();
            let mut positions_kthth_3_2 = Vec::new();
            let mut positions_kthu_3_3 = Vec::new();
            let mut positions_kthth_3_3 = Vec::new();
            let mut positions_kthu_3_4 = Vec::new();
            let mut positions_kthth_3_4 = Vec::new();

            let mut positions_kuu_4_1 = Vec::new();
            let mut positions_kuth_4_1 = Vec::new();
            let mut positions_kuu_4_2 = Vec::new();
            let mut positions_kuth_4_2 = Vec::new();
            let mut positions_kuu_4_3 = Vec::new();
            let mut positions_kuth_4_3 = Vec::new();
            let mut positions_kuu_4_4 = Vec::new();
            let mut positions_kuth_4_4 = Vec::new();

            let mut positions_kthu_4_1 = Vec::new();
            let mut positions_kthth_4_1 = Vec::new();
            let mut positions_kthu_4_2 = Vec::new();
            let mut positions_kthth_4_2 = Vec::new();
            let mut positions_kthu_4_3 = Vec::new();
            let mut positions_kthth_4_3 = Vec::new();
            let mut positions_kthu_4_4 = Vec::new();
            let mut positions_kthth_4_4 = Vec::new();

            let half = QuadFullPlateAuxFunctions::<T, V>::node_dof() / T::from(2u8);
            let one = QuadFullPlateAuxFunctions::<T, V>::node_dof();
            let one_and_half = QuadFullPlateAuxFunctions::<T, V>::node_dof() + 
                QuadFullPlateAuxFunctions::<T, V>::node_dof() / T::from(2u8);
            let two = QuadFullPlateAuxFunctions::<T, V>::node_dof() * T::from(2u8);
            let two_and_half = QuadFullPlateAuxFunctions::<T, V>::node_dof() * T::from(2u8) + 
                QuadFullPlateAuxFunctions::<T, V>::node_dof() / T::from(2u8);
            let three = QuadFullPlateAuxFunctions::<T, V>::node_dof() * T::from(3u8);
            let three_and_half = QuadFullPlateAuxFunctions::<T, V>::node_dof() * T::from(3u8) + 
                QuadFullPlateAuxFunctions::<T, V>::node_dof() / T::from(2u8);
    
            let mut i = T::from(0u8);
            while i < rows_number * columns_number
            {
                let position = MatrixElementPosition::create(
                    i / columns_number, i % columns_number);
    
                let row = i / columns_number;
                let column = i % columns_number;
    
                if row < half && column < half
                {
                    positions_kuu_1_1.push(position);
                }
                else if row < half && column >= half && column < one
                {
                    positions_kuth_1_1.push(position);
                }
                else if row < half && column >= one && column < one_and_half
                {
                    positions_kuu_1_2.push(position);
                }
                else if row < half && column >= one_and_half && column < two
                {
                    positions_kuth_1_2.push(position);
                }
                else if row < half && column >= two && column < two_and_half
                {
                    positions_kuu_1_3.push(position);
                }
                else if row < half && column >= two_and_half && column < three
                {
                    positions_kuth_1_3.push(position);
                }
                else if row < half && column >= three && column < three_and_half
                {
                    positions_kuu_1_4.push(position);
                }
                else if row < half && column >= three_and_half
                {
                    positions_kuth_1_4.push(position);
                }

                else if row >= half && row < one && column < half
                {
                    positions_kthu_1_1.push(position);
                }
                else if row >= half && row < one && column >= half && column < one
                {
                    positions_kthth_1_1.push(position);
                }
                else if row >= half && row < one && column >= one && column < one_and_half
                {
                    positions_kthu_1_2.push(position);
                }
                else if row >= half && row < one && column >= one_and_half && column < two
                {
                    positions_kthth_1_2.push(position);
                }
                else if row >= half && row < one && column >= two && column < two_and_half
                {
                    positions_kthu_1_3.push(position);
                }
                else if row >= half && row < one && column >= two_and_half && column < three
                {
                    positions_kthth_1_3.push(position);
                }
                else if row >= half && row < one && column >= three && column < three_and_half
                {
                    positions_kthu_1_4.push(position);
                }
                else if row >= half && row < one && column >= three_and_half
                {
                    positions_kthth_1_4.push(position);
                }

                else if row >= one && row < one_and_half && column < half
                {
                    positions_kuu_2_1.push(position);
                }
                else if row >= one && row < one_and_half && column >= half && column < one
                {
                    positions_kuth_2_1.push(position);
                }
                else if row >= one && row < one_and_half && column >= one && column < one_and_half
                {
                    positions_kuu_2_2.push(position);
                }
                else if row >= one && row < one_and_half && column >= one_and_half && column < two
                {
                    positions_kuth_2_2.push(position);
                }
                else if row >= one && row < one_and_half && column >= two && column < two_and_half
                {
                    positions_kuu_2_3.push(position);
                }
                else if row >= one && row < one_and_half && column >= two_and_half && column < three
                {
                    positions_kuth_2_3.push(position);
                }
                else if row >= one && row < one_and_half && column >= three && column < three_and_half
                {
                    positions_kuu_2_4.push(position);
                }
                else if row >= one && row < one_and_half && column >= three_and_half
                {
                    positions_kuth_2_4.push(position);
                }

                else if row >= one_and_half && row < two && column < half
                {
                    positions_kthu_2_1.push(position);
                }
                else if row >= one_and_half && row < two && column >= half && column < one
                {
                    positions_kthth_2_1.push(position);
                }
                else if row >= one_and_half && row < two && column >= one && column < one_and_half
                {
                    positions_kthu_2_2.push(position);
                }
                else if row >= one_and_half && row < two && column >= one_and_half && column < two
                {
                    positions_kthth_2_2.push(position);
                }
                else if row >= one_and_half && row < two && column >= two && column < two_and_half
                {
                    positions_kthu_2_3.push(position);
                }
                else if row >= one_and_half && row < two && column >= two_and_half && column < three
                {
                    positions_kthth_2_3.push(position);
                }
                else if row >= one_and_half && row < two && column >= three && column < three_and_half
                {
                    positions_kthu_2_4.push(position);
                }
                else if row >= one_and_half && row < two && column >= three_and_half
                {
                    positions_kthth_2_4.push(position);
                }

                else if row >= two && row < two_and_half && column < half
                {
                    positions_kuu_3_1.push(position);
                }
                else if row >= two && row < two_and_half && column >= half && column < one
                {
                    positions_kuth_3_1.push(position);
                }
                else if row >= two && row < two_and_half && column >= one && column < one_and_half
                {
                    positions_kuu_3_2.push(position);
                }
                else if row >= two && row < two_and_half && column >= one_and_half && column < two
                {
                    positions_kuth_3_2.push(position);
                }
                else if row >= two && row < two_and_half && column >= two && column < two_and_half
                {
                    positions_kuu_3_3.push(position);
                }
                else if row >= two && row < two_and_half && column >= two_and_half && column < three
                {
                    positions_kuth_3_3.push(position);
                }
                else if row >= two && row < two_and_half && column >= three && column < three_and_half
                {
                    positions_kuu_3_4.push(position);
                }
                else if row >= two && row < two_and_half && column >= three_and_half
                {
                    positions_kuth_3_4.push(position);
                }

                else if row >= two_and_half && row < three && column < half
                {
                    positions_kthu_3_1.push(position);
                }
                else if row >= two_and_half && row < three && column >= half && column < one
                {
                    positions_kthth_3_1.push(position);
                }
                else if row >= two_and_half && row < three && column >= one && column < one_and_half
                {
                    positions_kthu_3_2.push(position);
                }
                else if row >= two_and_half && row < three && column >= one_and_half && column < two
                {
                    positions_kthth_3_2.push(position);
                }
                else if row >= two_and_half && row < three && column >= two && column < two_and_half
                {
                    positions_kthu_3_3.push(position);
                }
                else if row >= two_and_half && row < three && column >= two_and_half && column < three
                {
                    positions_kthth_3_3.push(position);
                }
                else if row >= two_and_half && row < three && column >= three && column < three_and_half
                {
                    positions_kthu_3_4.push(position);
                }
                else if row >= two_and_half && row < three && column >= three_and_half
                {
                    positions_kthth_3_4.push(position);
                }

                else if row >= three && row < three_and_half && column < half
                {
                    positions_kuu_4_1.push(position);
                }
                else if row >= three && row < three_and_half && column >= half && column < one
                {
                    positions_kuth_4_1.push(position);
                }
                else if row >= three && row < three_and_half && column >= one && column < one_and_half
                {
                    positions_kuu_4_2.push(position);
                }
                else if row >= three && row < three_and_half && column >= one_and_half && column < two
                {
                    positions_kuth_4_2.push(position);
                }
                else if row >= three && row < three_and_half && column >= two && column < two_and_half
                {
                    positions_kuu_4_3.push(position);
                }
                else if row >= three && row < three_and_half && column >= two_and_half && column < three
                {
                    positions_kuth_4_3.push(position);
                }
                else if row >= three && row < three_and_half && column >= three && column < three_and_half
                {
                    positions_kuu_4_4.push(position);
                }
                else if row >= three && row < three_and_half && column >= three_and_half
                {
                    positions_kuth_4_4.push(position);
                }

                else if row >= three_and_half && column < half
                {
                    positions_kthu_4_1.push(position);
                }
                else if row >= three_and_half && column >= half && column < one
                {
                    positions_kthth_4_1.push(position);
                }
                else if row >= three_and_half && column >= one && column < one_and_half
                {
                    positions_kthu_4_2.push(position);
                }
                else if row >= three_and_half && column >= one_and_half && column < two
                {
                    positions_kthth_4_2.push(position);
                }
                else if row >= three_and_half && column >= two && column < two_and_half
                {
                    positions_kthu_4_3.push(position);
                }
                else if row >= three_and_half && column >= two_and_half && column < three
                {
                    positions_kthth_4_3.push(position);
                }
                else if row >= three_and_half && column >= three && column < three_and_half
                {
                    positions_kthu_4_4.push(position);
                }
                else 
                {
                    positions_kthth_4_4.push(position);
                }
                i += T::from(1u8);
            }
    
            vec![
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                    number_2: self.node_1_number, positions: positions_kuu_1_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_1_number,
                    number_2: self.node_1_number, positions: positions_kuth_1_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                    number_2: self.node_2_number, positions: positions_kuu_1_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_1_number,
                    number_2: self.node_2_number, positions: positions_kuth_1_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                    number_2: self.node_3_number, positions: positions_kuu_1_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_1_number,
                    number_2: self.node_3_number, positions: positions_kuth_1_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                    number_2: self.node_4_number, positions: positions_kuu_1_4, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_1_number,
                    number_2: self.node_4_number, positions: positions_kuth_1_4, },

                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_1_number,
                    number_2: self.node_1_number, positions: positions_kthu_1_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_1_number,
                    number_2: self.node_1_number, positions: positions_kthth_1_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_1_number,
                    number_2: self.node_2_number, positions: positions_kthu_1_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_1_number,
                    number_2: self.node_2_number, positions: positions_kthth_1_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_1_number,
                    number_2: self.node_3_number, positions: positions_kthu_1_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_1_number,
                    number_2: self.node_3_number, positions: positions_kthth_1_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_1_number,
                    number_2: self.node_4_number, positions: positions_kthu_1_4, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_1_number,
                    number_2: self.node_4_number, positions: positions_kthth_1_4, },
    
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                    number_2: self.node_1_number, positions: positions_kuu_2_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_2_number,
                    number_2: self.node_1_number, positions: positions_kuth_2_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                    number_2: self.node_2_number, positions: positions_kuu_2_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_2_number,
                    number_2: self.node_2_number, positions: positions_kuth_2_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                    number_2: self.node_3_number, positions: positions_kuu_2_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_2_number,
                    number_2: self.node_3_number, positions: positions_kuth_2_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                    number_2: self.node_4_number, positions: positions_kuu_2_4, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_2_number,
                    number_2: self.node_4_number, positions: positions_kuth_2_4, },

                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_2_number,
                    number_2: self.node_1_number, positions: positions_kthu_2_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_2_number,
                    number_2: self.node_1_number, positions: positions_kthth_2_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_2_number,
                    number_2: self.node_2_number, positions: positions_kthu_2_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_2_number,
                    number_2: self.node_2_number, positions: positions_kthth_2_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_2_number,
                    number_2: self.node_3_number, positions: positions_kthu_2_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_2_number,
                    number_2: self.node_3_number, positions: positions_kthth_2_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_2_number,
                    number_2: self.node_4_number, positions: positions_kthu_2_4, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_2_number,
                    number_2: self.node_4_number, positions: positions_kthth_2_4, },

                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_3_number,
                    number_2: self.node_1_number, positions: positions_kuu_3_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_3_number,
                    number_2: self.node_1_number, positions: positions_kuth_3_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_3_number,
                    number_2: self.node_2_number, positions: positions_kuu_3_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_3_number,
                    number_2: self.node_2_number, positions: positions_kuth_3_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_3_number,
                    number_2: self.node_3_number, positions: positions_kuu_3_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_3_number,
                    number_2: self.node_3_number, positions: positions_kuth_3_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_3_number,
                    number_2: self.node_4_number, positions: positions_kuu_3_4, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_3_number,
                    number_2: self.node_4_number, positions: positions_kuth_3_4, },

                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_3_number,
                    number_2: self.node_1_number, positions: positions_kthu_3_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_3_number,
                    number_2: self.node_1_number, positions: positions_kthth_3_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_3_number,
                    number_2: self.node_2_number, positions: positions_kthu_3_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_3_number,
                    number_2: self.node_2_number, positions: positions_kthth_3_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_3_number,
                    number_2: self.node_3_number, positions: positions_kthu_3_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_3_number,
                    number_2: self.node_3_number, positions: positions_kthth_3_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_3_number,
                    number_2: self.node_4_number, positions: positions_kthu_3_4, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_3_number,
                    number_2: self.node_4_number, positions: positions_kthth_3_4, },

                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_4_number,
                    number_2: self.node_1_number, positions: positions_kuu_4_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_4_number,
                    number_2: self.node_1_number, positions: positions_kuth_4_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_4_number,
                    number_2: self.node_2_number, positions: positions_kuu_4_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_4_number,
                    number_2: self.node_2_number, positions: positions_kuth_4_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_4_number,
                    number_2: self.node_3_number, positions: positions_kuu_4_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_4_number,
                    number_2: self.node_3_number, positions: positions_kuth_4_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_4_number,
                    number_2: self.node_4_number, positions: positions_kuu_4_4, },
                StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_4_number,
                    number_2: self.node_4_number, positions: positions_kuth_4_4, },

                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_4_number,
                    number_2: self.node_1_number, positions: positions_kthu_4_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_4_number,
                    number_2: self.node_1_number, positions: positions_kthth_4_1, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_4_number,
                    number_2: self.node_2_number, positions: positions_kthu_4_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_4_number,
                    number_2: self.node_2_number, positions: positions_kthth_4_2, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_4_number,
                    number_2: self.node_3_number, positions: positions_kthu_4_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_4_number,
                    number_2: self.node_3_number, positions: positions_kthth_4_3, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_4_number,
                    number_2: self.node_4_number, positions: positions_kthu_4_4, },
                StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_4_number,
                    number_2: self.node_4_number, positions: positions_kthth_4_4, },
            ]
    }


    fn is_node_belongs_to_element(&self, node_number: T) -> bool
    {
        self.node_1_number == node_number || self.node_2_number == node_number ||
        self.node_3_number == node_number || self.node_4_number == node_number
    }


    fn refresh(&mut self, tolerance: V, ref_nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        let rotation_matrix = QuadFullPlateAuxFunctions::rotation_matrix(
            self.node_2_number, self.node_3_number, self.node_4_number, tolerance, ref_nodes)?;

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(),
            QuadFullPlateAuxFunctions::<T, V>::nodes_number() * QuadFullPlateAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (PLATE4N4IP_NODES_NUMBER * PLATE_NODE_DOF).pow(2)], tolerance)?;

        for i in 0..self.state.integration_points.len()
        {
            let r = self.state.integration_points[i].r;
            let s = self.state.integration_points[i].s;
            let alpha = self.state.integration_points[i].weight_r * self.state.integration_points[i].weight_s;

            let matrix = QuadFullPlateAuxFunctions::local_stiffness_matrix(
                self.node_1_number, self.node_2_number, self.node_3_number, self.node_4_number, 
                self.young_modulus, self.poisson_ratio, self.thickness, self.shear_factor, alpha, r, s, ref_nodes,
                &rotation_matrix, tolerance)?;

            local_stiffness_matrix = local_stiffness_matrix.add_matrix(&matrix)
                .map_err(|e| format!("Plate4n4ip: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        self.state.rotation_matrix = rotation_matrix;
        self.state.local_stiffness_matrix = local_stiffness_matrix;
        Ok(())
    }


    fn is_nodes_numbers_same(&self, nodes_numbers: Vec<T>) -> bool
    {
        nodes_numbers.iter().all(|node_number| self.is_node_belongs_to_element(*node_number))
    }


    fn extract_element_analysis_data(&self, global_displacements: &Displacements<T, V>,
        tolerance: V, ref_nodes: &HashMap<T, FENode<V>>) -> Result<ElementAnalysisData<T, V>, String>
    {
        let element_local_displacements =
            self.extract_local_displacements(global_displacements, tolerance)?;

        let c_matrix_multiplier_mem = self.young_modulus / (V::from(1f32) - self.poisson_ratio.my_powi(2));
        let mut c_matrix_mem = ExtendedMatrix::create(
            T::from(3u8), T::from(3u8), 
            vec![
                V::from(1f32), self.poisson_ratio, V::from(0f32),
                self.poisson_ratio, V::from(1f32), V::from(0f32),
                V::from(0f32), V::from(0f32), (V::from(1f32) - self.poisson_ratio) / V::from(2f32),
            ], 
            tolerance)?;
        c_matrix_mem.multiply_by_number(c_matrix_multiplier_mem);

        let c_matrix_multiplier_bend = self.young_modulus * self.thickness.my_powi(3) / 
            (V::from(12f32) * (V::from(1f32) - self.poisson_ratio.my_powi(2)));
        let mut c_matrix_bend = ExtendedMatrix::create(
            T::from(3u8), T::from(3u8), 
            vec![
                V::from(1f32), self.poisson_ratio, V::from(0f32),
                self.poisson_ratio, V::from(1f32), V::from(0f32),
                V::from(0f32), V::from(0f32), (V::from(1f32) - self.poisson_ratio) / V::from(2f32),
            ], 
            tolerance)?;
        c_matrix_bend.multiply_by_number(c_matrix_multiplier_bend);

        let c_matrix_multiplier_shear = self.young_modulus * self.thickness * self.shear_factor / 
            (V::from(2f32) * (V::from(1f32) + self.poisson_ratio));
        let mut c_matrix_shear = ExtendedMatrix::create(
            T::from(2u8), T::from(2u8), 
            vec![
                V::from(1f32), V::from(0f32),
                V::from(0f32), V::from(1f32),
            ], 
            tolerance)?;
        c_matrix_shear.multiply_by_number(c_matrix_multiplier_shear);

        let mut forces_values = Vec::new();
        let mut forces_components = Vec::new();

        let mut force_x = V::from(0f32);
        let mut force_y = V::from(0f32);
        let mut force_xy = V::from(0f32);
        let mut moment_x = V::from(0f32);
        let mut moment_y = V::from(0f32);
        let mut moment_xy = V::from(0f32);
        let mut force_xz = V::from(0f32);
        let mut force_yz = V::from(0f32);

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

            let strain_displacement_matrix_mem_at_node = 
                QuadFullPlateAuxFunctions::strain_displacement_matrix_mem(self.node_1_number, self.node_2_number, 
                    self.node_3_number, self.node_4_number, r, s, ref_nodes, &self.state.rotation_matrix, tolerance)?;
            let strains_matrix_mem_at_node = 
                strain_displacement_matrix_mem_at_node.multiply_by_matrix(&element_local_displacements)?;

            let stresses_matrix_mem_at_node = c_matrix_mem.multiply_by_matrix(&strains_matrix_mem_at_node)?;
            let stresses_mem_at_node = 
                QuadFullPlateAuxFunctions::extract_column_matrix_values(&stresses_matrix_mem_at_node)?;
            for j in 0..3
            {
                if j == 0 { force_x += stresses_mem_at_node[j] * self.thickness; }
                if j == 1 { force_y += stresses_mem_at_node[j] * self.thickness; }
                if j == 2 { force_xy += stresses_mem_at_node[j] * self.thickness; }
            }

            let strain_displacement_matrix_bend_at_node = 
                QuadFullPlateAuxFunctions::strain_displacement_matrix_plate_bending(self.node_1_number, self.node_2_number, 
                    self.node_3_number, self.node_4_number, r, s, ref_nodes, &self.state.rotation_matrix, tolerance)?;
            let strains_matrix_bend_at_node = 
                strain_displacement_matrix_bend_at_node.multiply_by_matrix(&element_local_displacements)?;
            let stresses_matrix_bend_at_node = c_matrix_bend.multiply_by_matrix(&strains_matrix_bend_at_node)?;
            let stresses_bend_at_node = 
                QuadFullPlateAuxFunctions::extract_column_matrix_values(&stresses_matrix_bend_at_node)?;
            for k in 0..3
            {
                if k == 0 { moment_x += stresses_bend_at_node[k]; }
                if k == 1 { moment_y += stresses_bend_at_node[k]; }
                if k == 2 { moment_xy += stresses_bend_at_node[k]; }
            }

            let strain_displacement_matrix_shear_at_node = 
                QuadFullPlateAuxFunctions::strain_displacement_matrix_plate_shear(self.node_1_number, self.node_2_number, 
                    self.node_3_number, self.node_4_number, r, s, ref_nodes, &self.state.rotation_matrix, tolerance)?;
            let strains_matrix_shear_at_node = 
                strain_displacement_matrix_shear_at_node.multiply_by_matrix(&element_local_displacements)?;
            let stresses_matrix_shear_at_node = c_matrix_shear.multiply_by_matrix(&strains_matrix_shear_at_node)?;
            let stresses_shear_at_node = 
                QuadFullPlateAuxFunctions::extract_column_matrix_values(&stresses_matrix_shear_at_node)?;
            for m in 0..2
            {
                if m == 0 { force_xz += stresses_shear_at_node[m]; }
                if m == 1 { force_yz += stresses_shear_at_node[m]; }
            }
        }

        forces_values.push(force_x / V::from(4f32));
        forces_components.push(ForceComponent::MembraneForceX);
        forces_values.push(force_y / V::from(4f32));
        forces_components.push(ForceComponent::MembraneForceY);
        forces_values.push(force_xy / V::from(4f32));
        forces_components.push(ForceComponent::MembraneForceXY);
        forces_values.push(force_xz / V::from(4f32));
        forces_components.push(ForceComponent::ShearForceXZ);
        forces_values.push(force_yz / V::from(4f32));
        forces_components.push(ForceComponent::ShearForceYZ);
        forces_values.push(moment_x / V::from(4f32));
        forces_components.push(ForceComponent::MomentX);
        forces_values.push(moment_y / V::from(4f32));
        forces_components.push(ForceComponent::MomentY);
        forces_values.push(moment_xy / V::from(4f32));
        forces_components.push(ForceComponent::MomentXY);
    
        let element_forces = ElementForces::create(forces_values, forces_components);

        let element_analysis_data = ElementAnalysisData::create(
            None, None,
            Some(element_forces), None);
        Ok(element_analysis_data)
    }


    
    fn copy_nodes_numbers(&self) -> Vec<T>
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
}
