use std::hash::Hash;
use std::fmt::Debug;
use std::ops::{Sub, Mul, Add, Div, Rem, SubAssign, AddAssign, MulAssign};
use std::collections::HashMap;
use std::f32::consts::PI;

use extended_matrix::basic_matrix::basic_matrix::MatrixElementPosition;
use extended_matrix::extended_matrix::ExtendedMatrix;

use crate::fem::finite_elements::finite_element::{FiniteElementTrait, FEType};
use crate::fem::finite_elements::fe_node::FENode;

use crate::fem::global_analysis::fe_stiffness::{StiffnessGroup, StiffnessType};
use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData};
use crate::fem::global_analysis::fe_global_analysis_result::Displacements;

use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::element_analysis::fe_element_analysis_result::{ElementAnalysisData, ElementForces, NodalForces};

use crate::fem::finite_elements::beam::beam_aux_functions::BeamAuxFunctions;

use crate::my_float::MyFloatTrait;

use crate::fem::finite_elements::beam::consts::{ BEAM_NODE_DOF, BEAM2N1IPT_NODES_NUMBER };
use crate::fem::finite_elements::functions::extract_unique_elements_of_rotation_matrix;


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


pub struct Beam2n1ipT<T, V>
{
    node_1_number: T,
    node_2_number: T,
    young_modulus: V,
    poisson_ratio: V,
    area: V,
    i11: V,
    i22: V,
    angle: V,
    shear_factor: V,
    it: V,
    local_axis_1_direction: [V; 3],
    state: State<T, V>,
}


impl<T, V> Beam2n1ipT<T, V>
    where T: Copy + PartialOrd + Add<Output = T> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + AddAssign +
             From<u8> + 'static,
          V: Copy + Into<f64> + Sub<Output = V> + Mul<Output = V> + From<f32> + Add<Output = V> +
             Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign + SubAssign +
             MyFloatTrait + PartialOrd + MyFloatTrait<Other = V> + 'static
{
    pub fn create(node_1_number: T, node_2_number: T, young_modulus: V, poisson_ratio: V, area: V,
    i11_init: V, i22_init: V, i12_init: V, it: V, shear_factor: V,
    local_axis_1_direction: [V; 3], tolerance: V, nodes: &HashMap<T, FENode<V>>)
    -> Result<Self, String>
    {
        let mut angle = if i11_init == i22_init
            {
                V::from(0f32)
            }
            else
            {
                (V::from(2f32) * i12_init / (i22_init - i11_init)).my_atan() /
                V::from(2f32)
            };

        let mut i11 = i11_init * (angle.my_cos()).my_powi(2) +
            i22_init * (angle.my_sin()).my_powi(2) -
            i12_init * (V::from(2f32) * angle).my_sin();

        let mut i22 = i11_init * (angle.my_sin()).my_powi(2) +
            i22_init * (angle.my_cos()).my_powi(2) +
            i12_init * (V::from(2f32) * angle).my_sin();

        let mut i = 1;
        while i11 < i22
        {
            angle = ((V::from(2f32) * i12_init / (i22_init - i11_init)).my_atan() +
                V::from(PI) * V::from(i as f32)) / V::from(2f32);
            i11 = i11_init * (angle.my_cos()).my_powi(2) +
                i22_init * (angle.my_sin()).my_powi(2) -
                i12_init * (V::from(2f32) * angle).my_sin();
            i22 = i11_init * (angle.my_sin()).my_powi(2) +
                i22_init * (angle.my_cos()).my_powi(2) +
                i12_init * (V::from(2f32) * angle).my_sin();
            i += 1;
        }

        let integration_point_1 = IntegrationPoint {
            r: V::from(0f32), weight: V::from(2f32) };

        let rotation_matrix =
            BeamAuxFunctions::<T, V>::rotation_matrix(node_1_number, node_2_number,
                &local_axis_1_direction, angle, tolerance, nodes)?;

        let integration_points = vec![integration_point_1];

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (BEAM2N1IPT_NODES_NUMBER * BEAM_NODE_DOF).pow(2)], tolerance);

        for integration_point in &integration_points
        {
            let matrix = BeamAuxFunctions::local_stiffness_matrix(
                node_1_number, node_2_number, young_modulus, poisson_ratio, area, shear_factor,
                it, i11, i22, integration_point.weight, integration_point.r,
                tolerance, nodes)?;

            local_stiffness_matrix = local_stiffness_matrix.add_matrix(&matrix)
                .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        let mut nodes_dof_parameters =
            BeamAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;

        let node_2_dof_parameters =
            BeamAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;

        nodes_dof_parameters.extend(node_2_dof_parameters);

        let state = State::create(rotation_matrix, integration_points,
            local_stiffness_matrix, nodes_dof_parameters);

        Ok(Beam2n1ipT { node_1_number, node_2_number, young_modulus, poisson_ratio,
            area, i11, i22, angle, shear_factor, it, local_axis_1_direction, state })
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


impl<T, V> FiniteElementTrait<T, V> for Beam2n1ipT<T, V>
    where T: Copy + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Rem<Output = T> +
             Mul<Output = T> + Eq + Hash + Debug + SubAssign + PartialOrd + AddAssign +
             From<u8> + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             Into<f64> + SubAssign + AddAssign + MulAssign + PartialEq + Debug +
             MyFloatTrait + PartialOrd + From<f32> + MyFloatTrait<Other = V> + 'static,
{
    fn update(&mut self, nodes_numbers: Vec<T>, properties: Vec<V>, tolerance: V,
        nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        let node_1_number = nodes_numbers[0];

        let node_2_number = nodes_numbers[1];

        let young_modulus = properties[0];

        let poisson_ratio = properties[1];

        let area = properties[2];

        let i11_init = properties[3];

        let i22_init = properties[4];

        let i12_init = properties[5];

        let mut angle = if i11_init == i22_init
            {
                V::from(0f32)
            }
            else
            {
                (V::from(2f32) * i12_init / (i22_init - i11_init)).my_atan() /
                V::from(2f32)
            };

        let mut i11 = i11_init * (angle.my_cos()).my_powi(2) +
            i22_init * (angle.my_sin()).my_powi(2) -
            i12_init * (V::from(2f32) * angle).my_sin();

        let mut i22 = i11_init * (angle.my_sin()).my_powi(2) +
            i22_init * (angle.my_cos()).my_powi(2) +
            i12_init * (V::from(2f32) * angle).my_sin();

        let mut i = 1;
        while i11 < i22
        {
            angle = ((V::from(2f32) * i12_init / (i22_init - i11_init)).my_atan() +
                V::from(PI) * V::from(i as f32)) / V::from(2f32);
            i11 = i11_init * (angle.my_cos()).my_powi(2) +
                i22_init * (angle.my_sin()).my_powi(2) -
                i12_init * (V::from(2f32) * angle).my_sin();
            i22 = i11_init * (angle.my_sin()).my_powi(2) +
                i22_init * (angle.my_cos()).my_powi(2) +
                i12_init * (V::from(2f32) * angle).my_sin();
            i += 1;
        }

        let it = properties[6];

        let shear_factor = properties[7];

        let local_axis_1_direction = [properties[8], properties[9], properties[10]];

        let rotation_matrix =
            BeamAuxFunctions::rotation_matrix(node_1_number, node_2_number,
                &local_axis_1_direction, angle, tolerance, nodes)?;

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (BEAM2N1IPT_NODES_NUMBER * BEAM_NODE_DOF).pow(2)], tolerance);

        for integration_point in &self.state.integration_points
        {
            let matrix = BeamAuxFunctions::local_stiffness_matrix(
                node_1_number, node_2_number, young_modulus, poisson_ratio, area, shear_factor,
                it, i11, i22, integration_point.weight, integration_point.r,
                tolerance, nodes)?;
            local_stiffness_matrix = local_stiffness_matrix.add_matrix(&matrix)
                .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
        }

        let mut nodes_dof_parameters =
            BeamAuxFunctions::<T, V>::compose_node_dof_parameters(node_1_number)?;

        let node_2_dof_parameters =
            BeamAuxFunctions::<T, V>::compose_node_dof_parameters(node_2_number)?;

        nodes_dof_parameters.extend(node_2_dof_parameters);

        self.node_1_number = node_1_number;
        self.node_2_number = node_2_number;
        self.young_modulus = young_modulus;
        self.poisson_ratio = poisson_ratio;
        self.area = area;
        self.i11 = i11;
        self.i22 = i22;
        self.angle = angle;
        self.shear_factor = shear_factor;
        self.it = it;
        self.local_axis_1_direction = local_axis_1_direction;
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
        Err("Beam2n2ipT: Stiffness matrix cannot be extracted!")
    }


    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>
    {
        let (rows_number, columns_number) =
            (BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
             BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof());

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

        let mut i = T::from(0u8);
        while i < rows_number * columns_number
        {
            let position = MatrixElementPosition::create(
                i / columns_number, i % columns_number);

            let row = i / columns_number;
            let column = i % columns_number;

            if row < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)
            {
                positions_kuu_1_1.push(position);
            }
            else if row < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof()
            {
                positions_kuth_1_1.push(position);
            }
            else if row < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() &&
                column < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kuu_1_2.push(position);
            }
            else if row < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kuth_1_2.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                row < BeamAuxFunctions::<T, V>::node_dof() &&
                column < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)
            {
                positions_kthu_1_1.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                row < BeamAuxFunctions::<T, V>::node_dof() &&
                column >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof()
            {
                positions_kthth_1_1.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                row < BeamAuxFunctions::<T, V>::node_dof() &&
                column >= BeamAuxFunctions::<T, V>::node_dof() &&
                column < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kthu_1_2.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                row < BeamAuxFunctions::<T, V>::node_dof() &&
                column >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kthth_1_2.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() &&
                row < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)
            {
                positions_kuu_2_1.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() &&
                row < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof()
            {
                positions_kuth_2_1.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() &&
                row < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() &&
                column < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kuu_2_2.push(position);
            }
            else if row >= BeamAuxFunctions::<T, V>::node_dof() &&
                row < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kuth_2_2.push(position);
            }
            else if row >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column < BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)
            {
                positions_kthu_2_1.push(position);
            }
            else if row >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8) &&
                column < BeamAuxFunctions::<T, V>::node_dof()
            {
                positions_kthth_2_1.push(position);
            }
            else if row >= (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8)) &&
                column >= BeamAuxFunctions::<T, V>::node_dof() &&
                column < (BeamAuxFunctions::<T, V>::node_dof() +
                    BeamAuxFunctions::<T, V>::node_dof() / T::from(2u8))
            {
                positions_kthu_2_2.push(position);
            }
            else
            {
                positions_kthth_2_2.push(position);
            }
            i += T::from(1u8);
        }

        vec![StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                number_2: self.node_1_number, positions: positions_kuu_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_1_number,
                number_2: self.node_1_number, positions: positions_kuth_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_1_number,
                number_2: self.node_2_number, positions: positions_kuu_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_1_number,
                number_2: self.node_2_number, positions: positions_kuth_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_1_number,
                number_2: self.node_1_number, positions: positions_kthu_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_1_number,
                number_2: self.node_1_number, positions: positions_kthth_1_1, },
             StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_1_number,
                number_2: self.node_2_number, positions: positions_kthu_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_1_number,
                number_2: self.node_2_number, positions: positions_kthth_1_2, },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                number_2: self.node_1_number, positions: positions_kuu_2_1 },
             StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_2_number,
                number_2: self.node_1_number, positions: positions_kuth_2_1 },
             StiffnessGroup { stiffness_type: StiffnessType::Kuu, number_1: self.node_2_number,
                number_2: self.node_2_number, positions: positions_kuu_2_2 },
             StiffnessGroup { stiffness_type: StiffnessType::Kuth, number_1: self.node_2_number,
                number_2: self.node_2_number, positions: positions_kuth_2_2 },
             StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_2_number,
                number_2: self.node_1_number, positions: positions_kthu_2_1 },
             StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_2_number,
                number_2: self.node_1_number, positions: positions_kthth_2_1 },
             StiffnessGroup { stiffness_type: StiffnessType::Kthu, number_1: self.node_2_number,
                number_2: self.node_2_number, positions: positions_kthu_2_2 },
             StiffnessGroup { stiffness_type: StiffnessType::Kthth, number_1: self.node_2_number,
                number_2: self.node_2_number, positions: positions_kthth_2_2 },
        ]
    }


    fn is_node_belongs_to_element(&self, node_number: T) -> bool
    {
        self.node_1_number == node_number || self.node_2_number == node_number
    }


    fn refresh(&mut self, tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<(), String>
    {
        let rotation_matrix =
            BeamAuxFunctions::rotation_matrix(self.node_1_number, self.node_2_number,
                &self.local_axis_1_direction, self.angle, tolerance, nodes)?;

        let mut local_stiffness_matrix = ExtendedMatrix::create(
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            BeamAuxFunctions::<T, V>::nodes_number() * BeamAuxFunctions::<T, V>::node_dof(),
            vec![V::from(0f32); (BEAM2N1IPT_NODES_NUMBER * BEAM_NODE_DOF).pow(2)], tolerance);

        for integration_point in self.state.integration_points.iter()
        {
            let matrix = BeamAuxFunctions::local_stiffness_matrix(
                self.node_1_number, self.node_2_number, self.young_modulus, self.poisson_ratio,
                self.area, self.shear_factor, self.it, self.i11, self.i22,
                integration_point.weight, integration_point.r, tolerance, nodes)?;
            local_stiffness_matrix = local_stiffness_matrix.add_matrix(&matrix)
                .map_err(|e| format!("Beam2n2ipT: Local stiffness matrix could not be \
                    calculated! Reason: {}", e))?;
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


    fn extract_element_analysis_data(&self, global_displacements: &Displacements<T, V>,
        tolerance: V, nodes: &HashMap<T, FENode<V>>) -> Result<ElementAnalysisData<T, V>, String>
    {
        let element_local_displacements =
            self.extract_local_displacements(global_displacements, tolerance)?;

        let r = self.state.integration_points[0].r;
        let mut forces_values = Vec::new();
        let mut forces_components = Vec::new();

        let strain_displacement_matrix_u =
            BeamAuxFunctions::strain_displacement_matrix_u(
                self.node_1_number, self.node_2_number, r, tolerance, nodes);
        let strains_matrix_u =
            strain_displacement_matrix_u.multiply_by_matrix(&element_local_displacements)?;
        let force_x = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_u)[0] * self.area * self.young_modulus;
        forces_components.push(ForceComponent::ForceX);
        forces_values.push(force_x);

        let strain_displacement_matrix_v =
            BeamAuxFunctions::strain_displacement_matrix_v(
                self.node_1_number, self.node_2_number, r, tolerance, nodes)?;
        let strains_matrix_v =
            strain_displacement_matrix_v.multiply_by_matrix(&element_local_displacements)?;
        let shear_modulus = self.young_modulus /
            (V::from(2f32) * (V::from(1f32) + self.poisson_ratio));
        let force_y = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_v)[0] * shear_modulus * self.area * self.shear_factor;
        forces_components.push(ForceComponent::ForceY);
        forces_values.push(force_y);

        let strain_displacement_matrix_w =
            BeamAuxFunctions::strain_displacement_matrix_w(
                self.node_1_number, self.node_2_number, r, tolerance, nodes)?;
        let strains_matrix_w =
            strain_displacement_matrix_w.multiply_by_matrix(&element_local_displacements)?;
        let force_z = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_w)[0] * shear_modulus * self.area * self.shear_factor;
        forces_components.push(ForceComponent::ForceZ);
        forces_values.push(force_z);

        let strain_displacement_matrix_thu =
            BeamAuxFunctions::strain_displacement_matrix_thu(
                self.node_1_number, self.node_2_number, r, tolerance, nodes);
        let strains_matrix_thu =
            strain_displacement_matrix_thu.multiply_by_matrix(&element_local_displacements)?;
        let moment_x = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_thu)[0] * shear_modulus * self.it;
        forces_components.push(ForceComponent::MomentX);
        forces_values.push(moment_x);

        let length = BeamAuxFunctions::length(self.node_1_number, self.node_2_number, nodes);

        let mut forces_values_for_node_1 = Vec::new();
        let mut forces_components_for_node_1 = Vec::new();
        let mut forces_values_for_node_2 = Vec::new();
        let mut forces_components_for_node_2 = Vec::new();

        let strain_displacement_matrix_thv =
            BeamAuxFunctions::strain_displacement_matrix_thv(
                self.node_1_number, self.node_2_number, r, tolerance, nodes);
        let strains_matrix_thv =
            strain_displacement_matrix_thv.multiply_by_matrix(&element_local_displacements)?;
        let moment_y_average = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_thv)[0] * self.young_modulus * self.i22;
        forces_components.push(ForceComponent::MomentY);
        forces_values.push(moment_y_average);
        let moment_y_at_node_1 = moment_y_average + length * force_z / V::from(2f32);
        let moment_y_at_node_2 = moment_y_average - length * force_z / V::from(2f32);
        forces_components_for_node_1.push(ForceComponent::MomentY);
        forces_components_for_node_2.push(ForceComponent::MomentY);
        forces_values_for_node_1.push(moment_y_at_node_1);
        forces_values_for_node_2.push(moment_y_at_node_2);

        let strain_displacement_matrix_thw =
            BeamAuxFunctions::strain_displacement_matrix_thw(
                self.node_1_number, self.node_2_number, r, tolerance, nodes);
        let strains_matrix_thw =
            strain_displacement_matrix_thw.multiply_by_matrix(&element_local_displacements)?;
        let moment_z_average = BeamAuxFunctions::extract_column_matrix_values(
            &strains_matrix_thw)[0] * self.young_modulus * self.i11;
        forces_components.push(ForceComponent::MomentZ);
        forces_values.push(moment_z_average);
        let moment_z_at_node_1 = moment_z_average + length * force_y / V::from(2f32);
        let moment_z_at_node_2 = moment_z_average - length * force_y / V::from(2f32);
        forces_components_for_node_1.push(ForceComponent::MomentZ);
        forces_components_for_node_2.push(ForceComponent::MomentZ);
        forces_values_for_node_1.push(moment_z_at_node_1);
        forces_values_for_node_2.push(moment_z_at_node_2);

        let element_forces = ElementForces::create(forces_values,
            forces_components);

        let nodal_forces_for_node_1 = NodalForces::create(
            forces_values_for_node_1, forces_components_for_node_1);
        let nodal_forces_for_node_2 = NodalForces::create(
            forces_values_for_node_2, forces_components_for_node_2);

        let mut nodal_forces = HashMap::new();
        nodal_forces.insert(self.node_1_number, nodal_forces_for_node_1);
        nodal_forces.insert(self.node_2_number, nodal_forces_for_node_2);

        let element_analysis_data = ElementAnalysisData::create(
            None, None, Some(element_forces), Some(nodal_forces));
        Ok(element_analysis_data)

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


    fn extract_unique_elements_of_rotation_matrix(&self) -> Vec<V>
    {
        extract_unique_elements_of_rotation_matrix(&self.state.rotation_matrix)
    }
}
