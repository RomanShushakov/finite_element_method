use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign, DivAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::collections::{HashSet, HashMap};
use std::iter::FromIterator;

use extended_matrix::basic_matrix::basic_matrix::{MatrixElementPosition, ZerosRowColumn};
use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::extended_matrix::Operation;
use extended_matrix::functions::{extract_element_value, conversion_uint_into_usize};

use crate::fem::finite_elements::fe_node::FENode;
use crate::fem::finite_elements::finite_element::{FiniteElement, FEType};
use crate::fem::global_analysis::fe_stiffness::StiffnessGroup;
use crate::fem::global_analysis::fe_boundary_condition::{BoundaryCondition, BCType};
use crate::fem::global_analysis::fe_global_analysis_result::{GlobalAnalysisResult, Displacements};
use crate::fem::global_analysis::fe_dof_parameter_data::
{
    global_dof, DOFParameterData, GLOBAL_DOF, GlobalDOFParameter
};

use crate::fem::element_analysis::fe_element_analysis_result::{ElementAnalysisData, ElementsAnalysisResult};
use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
use crate::fem::element_analysis::beam::beam_element_nodal_forces::{NodalForces, BeamElementNodalForces};

use crate::fem::functions::{separate, compose_stiffness_sub_groups};

use crate::my_float::MyFloatTrait;


struct State<T, V>
{
    stiffness_groups: Vec<StiffnessGroup<T>>,
    nodes_dof_parameters_global: Vec<DOFParameterData<T>>,
    tolerance: V,
}


impl<T, V> State<T, V>
{
    fn create(stiffness_groups: Vec<StiffnessGroup<T>>,
        nodes_dof_parameters_global: Vec<DOFParameterData<T>>, tolerance: V,) -> Self
    {
        State { stiffness_groups, nodes_dof_parameters_global, tolerance }
    }
}


pub struct FEModel<T, V>
{
    pub nodes: HashMap<T, FENode<V>>,                           // Hashmap { node_number: Node }
    pub elements: HashMap<T, FiniteElement<T, V>>,              // Hashmap { element_number: FiniteElement }
    pub boundary_conditions: Vec<BoundaryCondition<T, V>>,
    state: State<T, V>,
}


impl<T, V> FEModel<T, V>
    where T: Copy + PartialEq + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash +
             SubAssign + Debug + Mul<Output = T> + PartialOrd + Add<Output = T> + From<u8> +
             AddAssign + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             PartialEq + Debug + AddAssign + MulAssign + SubAssign + Into<f64> + PartialOrd +
             MyFloatTrait + From<f32> + MyFloatTrait<Other = V> + DivAssign + 'static,
{
    pub fn create(tolerance: V) -> Self
    {
        let state = State::create(Vec::new(),
            Vec::new(), tolerance);

        FEModel { nodes: HashMap::new(), elements: HashMap::new(),
            boundary_conditions: Vec::new(), state }
    }


    fn update_stiffness_groups(&mut self) -> Result<(), &str>
    {
        let mut stiffness_groups = Vec::new();
        if self.nodes.len() < 2
        {
            self.state.stiffness_groups = Vec::new();
        }
        else
        {
            let mut nodes_numbers = Vec::new();
            for node_number in self.nodes.keys()
            {
                nodes_numbers.push(*node_number);
            }
            let mut position = T::from(0u8);

            let mut columns_number = T::from(0u8);
            (0..nodes_numbers.len()).for_each(|_| columns_number += T::from(1u8));

            for i in 1..nodes_numbers.len()
            {
                let mut v_lhs = nodes_numbers[0..i - 1].to_vec();
                let v_rhs = &nodes_numbers[i..];
                v_lhs.extend(v_rhs);
                let excluded = nodes_numbers[i - 1];
                for j in 0..v_lhs.len()
                {
                    if j + 1 == i
                    {
                        let stiffness_sub_groups =
                             compose_stiffness_sub_groups(position,
                            columns_number, nodes_numbers[j],
                            nodes_numbers[j])?;
                        stiffness_groups.extend(stiffness_sub_groups);
                        position += T::from(1u8);
                    }
                    let stiffness_sub_groups =
                         compose_stiffness_sub_groups(position,
                        columns_number, excluded,
                        v_lhs[j])?;
                    stiffness_groups.extend(stiffness_sub_groups);
                    position += T::from(1u8);
                }
            }

            for i in 0..nodes_numbers.len() - 1
            {
                let stiffness_sub_groups =
                     compose_stiffness_sub_groups(position,
                    columns_number,
                    nodes_numbers[nodes_numbers.len() - 1],
                    nodes_numbers[i])?;
                stiffness_groups.extend(stiffness_sub_groups);
                position += T::from(1u8);
            }
            let stiffness_sub_groups =
                 compose_stiffness_sub_groups(position,
                columns_number,
                nodes_numbers[nodes_numbers.len() - 1],
                nodes_numbers[nodes_numbers.len() - 1])?;
            stiffness_groups.extend(stiffness_sub_groups);
        }
        self.state.stiffness_groups = stiffness_groups;
        Ok(())
    }


    fn update_nodes_dof_parameters_global(&mut self) -> Result<(), &str>
    {
        let mut nodes_dof_parameters = Vec::new();
        for node_number in self.nodes.keys()
        {
            for dof in 0..GLOBAL_DOF
            {
                let dof_parameter =
                    GlobalDOFParameter::iterator().nth(dof)
                        .ok_or("FEModel: Could not find dof parameter!")?;
                let dof_parameter_data = DOFParameterData::create(*node_number,
                    *dof_parameter);
                nodes_dof_parameters.push(dof_parameter_data);
            }
        }
        self.state.nodes_dof_parameters_global = nodes_dof_parameters;
        Ok(())
    }


    pub fn add_node(&mut self, number: T, x: V, y: V, z: V) -> Result<(), String>
    {
        if self.nodes.contains_key(&number)
        {
            return Err(format!("FEModel: Node {:?} could not be added because node with same number \
                does already exist!", number));
        }

        if self.nodes.values().position(|node| node.is_coordinates_same(x, y, z)).is_some()
        {
            return Err(format!("FEModel: Node {:?} could not be added because node with the same \
                coordinates does already exist!", number));
        }

        let node = FENode::create(x, y, z);
        self.nodes.insert(number, node);
        self.update_stiffness_groups()?;
        Ok(())

    }


    pub fn update_node(&mut self, number: T, x: V, y: V, z: V) -> Result<(), String>
    {
        if self.nodes.iter().position(|(node_number, node)|
            *node_number != number && node.is_coordinates_same(x, y, z)).is_some()
        {
            return Err(format!("FEModel: Node {:?} could not be updated because the node with the \
                same coordinates does already exist!", number))
        }

        if let Some(node) = self.nodes.get_mut(&number)
        {
            node.update(x, y, z);
            for element in self.elements.values_mut()
                .filter(|element| element.is_node_belong_element(number))
            {
                element.refresh(self.state.tolerance, &self.nodes)?;
            }
            Ok(())
        }
        else
        {
            return Err(format!("FEModel: Node {:?} could not be updated because it does not \
                exist!", number));
        }
    }


    pub fn delete_node(&mut self, number: T) -> Result<(), String>
    {
        if !self.nodes.contains_key(&number)
        {
            return Err(format!("FEModel: Node {:?} could not be deleted because it does not exist!",
                number));
        }

        let mut elements_numbers_for_delete = Vec::new();
        for (element_number, element) in self.elements.iter()
        {
            if element.is_node_belong_element(number)
            {
                elements_numbers_for_delete.push(*element_number);
            }
        }

        for element_number in elements_numbers_for_delete
        {
            let _ = self.elements.remove(&element_number);
        }

        let _ = self.nodes.remove(&number);

        self.update_stiffness_groups()?;
        Ok(())
    }


    pub fn add_element(&mut self, element_number: T, element_type: FEType, nodes_numbers: Vec<T>,
        properties: Vec<V>) -> Result<(), String>
    {
        for (i, value) in properties.iter().enumerate()
        {
            if element_type != FEType::Beam2n1ipT
            {
                if *value <= V::from(0f32)
                {
                    return Err(format!("FEData: All properties values for element {:?} should be \
                        greater than zero!", element_number));
                }
            }
            else
            {
                if *value <= V::from(0f32) && [i != 5, i < 8].iter().all(|condition| *condition == true)
                {
                    return Err(format!("FEData: All properties values for element {:?} should be \
                        greater than zero!", element_number));
                }
            }
        }

        if self.elements.contains_key(&element_number)
        {
            return Err(format!("FEModel: Element {:?} could not be added! The element with the same \
             number does already exist!", element_number));
        }

        let nodes_numbers_set = HashSet::<T>::from_iter(
            nodes_numbers.iter().cloned());
        if nodes_numbers.len() != nodes_numbers_set.len()
        {
            return Err(format!("FEModel: Element {:?} could not be added! All nodes numbers \
                should be unique!", element_number));
        }

        if self.elements.values().position(|element|
            element.is_type_same(&element_type) &&
            element.is_nodes_numbers_same(nodes_numbers.clone())).is_some()
        {
            return Err(format!("FEModel: Element {:?} could not be added! The element with the same \
                type and with same nodes numbers does already exist!", element_number));
        }

        for node_number in nodes_numbers.iter()
        {
            if !self.nodes.contains_key(node_number)
            {
                return Err(format!("FEModel: Element {:?} could not be added! Node {:?} does not \
                    exist!", element_number, node_number));
            }
        }

        let element = FiniteElement::create(element_type, nodes_numbers,
            properties, self.state.tolerance, &self.nodes)?;
        self.elements.insert(element_number, element);

        Ok(())
    }


    pub fn update_element(&mut self, element_number: T, nodes_numbers: Vec<T>,
        properties: Vec<V>) -> Result<(), String>
    {
        for value in properties.iter()
        {
            if *value <= V::from(0f32)
            {
                return Err(format!("FEData: All properties values for element {:?} should be \
                    greater than zero!", element_number));
            }
        }

        let nodes_numbers_set = HashSet::<T>::from_iter(
            nodes_numbers.iter().cloned());
        if nodes_numbers.len() != nodes_numbers_set.len()
        {
            return Err(format!("FEModel: Element {:?} could not be updated! All nodes numbers \
                should be unique!", element_number));
        }

        for node_number in nodes_numbers.iter()
        {
            if !self.nodes.contains_key(node_number)
            {
                return Err(format!("FEModel: Element {:?} could not be updated! Node {:?} does not \
                    exist!", element_number, node_number));
            }
        }

        if let Some(element) = self.elements.get_mut(&element_number)
        {
            if element.is_nodes_numbers_same(nodes_numbers.clone())
            {
                return Err(format!("FEModel: Element {:?} could not be added! The element with \
                    the same nodes numbers does already exist!", element_number));
            }
            else
            {
                element.update(nodes_numbers, properties, self.state.tolerance, &self.nodes)?;
                Ok(())
            }
        }
        else
        {
            return Err(format!("FEModel: Element {:?} could not be updated because it does not \
                exist!", element_number));
        }
    }


    pub fn delete_element(&mut self, number: T) -> Result<(), String>
    {
        if self.elements.remove(&number).is_none()
        {
            return Err(format!("FEModel: Element {:?} could not be deleted because it does not \
                exist!", number))
        }
        Ok(())
    }


    fn compose_global_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
    {
        if self.elements.is_empty()
        {
            return Err("FEModel: Global stiffness matrix could not be composed because there are \
                no elements in the model!");
        }

        if self.nodes.keys().any(|node_number|
            self.elements.values().position(|element|
                element.is_node_belong_element(*node_number)).is_none())
        {
            return Err("FEModel: Global stiffness matrix could not be composed because there are \
                free nodes exist!");
        }

        let mut nodes_len_value = T::from(0u8);
        (0..self.nodes.len()).for_each(|_| nodes_len_value += T::from(1u8));

        let mut global_stiffness_matrix = ExtendedMatrix::<T, V>::create(
            nodes_len_value * global_dof::<T>(),
            nodes_len_value * global_dof::<T>(),
            vec![V::from(0f32); (self.nodes.len() * GLOBAL_DOF).pow(2)],
            self.state.tolerance);

        for element in self.elements.values()
        {
            let element_stiffness_matrix = element.extract_stiffness_matrix()?;

            let element_stiffness_groups = element.extract_stiffness_groups();
            for element_stiffness_group in element_stiffness_groups
            {
                if let Some(position) = self.state.stiffness_groups
                    .iter()
                    .position(|group|
                        { group.stiffness_type == element_stiffness_group.stiffness_type &&
                        group.number_1 == element_stiffness_group.number_1 &&
                        group.number_2 == element_stiffness_group.number_2 })
                {
                    global_stiffness_matrix.add_sub_matrix(
                        &element_stiffness_matrix,
                        &self.state.stiffness_groups[position].positions,
                        &element_stiffness_group.positions, self.state.tolerance);
                }
            }
        }
        Ok(global_stiffness_matrix)
    }


    pub fn add_bc(&mut self, bc_type: BCType, number: T, node_number: T,
        dof_parameter: GlobalDOFParameter, value: V) -> Result<(), String>
    {
        if self.boundary_conditions.iter().position(|bc|
            bc.is_number_same(number) && bc.is_type_same(bc_type)).is_some()
        {
            return Err(format!("FEModel: {} could not be added because the same {} number does \
                already exist!", bc_type.as_str(), bc_type.as_str().to_lowercase()));
        }

        if self.boundary_conditions.iter().position(|bc|
            bc.is_dof_parameter_data_same(dof_parameter, node_number)).is_some()
        {
            return Err(format!("FEModel: {} could not be added because the the force or \
                displacement with the same dof parameter data does already exist!",
                bc_type.as_str()));
        }

        if !self.nodes.contains_key(&node_number)
        {
            return Err(format!("FEModel: {} could not be added because the current node number \
                does not exist!", bc_type.as_str()));
        }

        let bc = BoundaryCondition::create(bc_type, number, node_number,
            dof_parameter, value);

        self.boundary_conditions.push(bc);
        Ok(())
    }


    pub fn update_bc(&mut self, bc_type: BCType, number: T, node_number: T,
        dof_parameter: GlobalDOFParameter, value: V) -> Result<(), String>
    {
        if !self.nodes.contains_key(&node_number)
        {
            return Err(format!("FEModel: {} could not be updated because the current node number \
                does not exist!", bc_type.as_str()));
        }

        if self.boundary_conditions.iter().position(|bc|
            (bc.is_dof_parameter_data_same(dof_parameter, node_number) &&
            !bc.is_number_same(number)) ||
            (bc.is_dof_parameter_data_same(dof_parameter, node_number) &&
            bc.is_number_same(number) && !bc.is_type_same(bc_type))).is_some()
        {
            return Err(format!("FEModel: {} could not be updated because the the force or \
                displacement with the same dof parameter data does already exist!",
                bc_type.as_str()));
        }

        if let Some(position) =  self.boundary_conditions.iter().position(|bc|
            bc.is_number_same(number) && bc.is_type_same(bc_type))
        {
            self.boundary_conditions[position].update(node_number, dof_parameter, value);
            Ok(())
        }
        else
        {
            Err(format!("FEModel: {} could not be updated because current {} number does not \
                exist!", bc_type.as_str(), bc_type.as_str().to_lowercase()))
        }
    }


    pub fn delete_bc(&mut self, bc_type: BCType, number: T) -> Result<(), String>
    {
        if let Some(position) =  self.boundary_conditions.iter().position(|bc|
            bc.is_number_same(number) && bc.is_type_same(bc_type))
        {
            self.boundary_conditions.remove(position);
            Ok(())
        }
        else
        {
            Err(format!("FEModel: {} could not be deleted because current {} number does not \
                exist!", bc_type.as_str(), bc_type.as_str().to_lowercase()))
        }
    }


    fn shrink_of_nodes_dof_parameters(&mut self, zeros_rows_columns: &Vec<ZerosRowColumn<T>>)
        -> Result<(), String>
    {
        for row_column in zeros_rows_columns
        {
            let mut row_column_as_index = 0usize;
            let mut n = T::from(0u8);
            while n < row_column.column()
            {
                row_column_as_index += 1usize;
                n += T::from(1u8);
            }

            let dof_parameter_data =
                self.state.nodes_dof_parameters_global.remove(row_column_as_index);

            if let Some(position) = self.boundary_conditions
                .iter()
                .position(|bc|
                    bc.is_dof_parameter_data_same(
                        dof_parameter_data.dof_parameter(),
                        dof_parameter_data.node_number()))
            {
                let bc_type = self.boundary_conditions[position].extract_bc_type();
                let dof_parameter = dof_parameter_data.dof_parameter();
                let node_number = dof_parameter_data.node_number();
                return Err(format!("FEModel: Model could not be analyzed because where are \
                    no stiffness to withstand {}::{:?} applied at node {:?}!", bc_type.as_str(),
                    dof_parameter, node_number))
            }
        }
        Ok(())
    }


    fn compose_separation_positions(&self, ub_rb_rows_numbers: &mut Vec<T>,
        separation_positions: &mut Vec<MatrixElementPosition<T>>)
    {
        for bc in &self.boundary_conditions
        {
            if bc.is_type_same(BCType::Displacement)
            {
                let mut row = T::from(0u8);
                for dof_parameter_data in
                    &self.state.nodes_dof_parameters_global
                {
                    if bc.is_dof_parameter_data_same(
                        dof_parameter_data.dof_parameter(),
                        dof_parameter_data.node_number())
                    {
                        separation_positions.push(
                            MatrixElementPosition::create(row, row));
                        ub_rb_rows_numbers.push(row);
                    }
                    row += T::from(1u8);
                }
            }
        }
    }


    fn compose_ua_ra_rows_numbers(&self, ub_rb_rows_numbers: &Vec<T>,
        ua_ra_rows_numbers: &mut Vec<T>)
    {
        let mut i = T::from(0u8);
        (0..self.state.nodes_dof_parameters_global.len()).for_each(|_|
            {
                if ub_rb_rows_numbers.iter().position(|n| *n == i).is_none()
                {
                    ua_ra_rows_numbers.push(i);
                }
                i += T::from(1u8);
            });
    }


    fn compose_matrix_by_rows_numbers(&self, rows_numbers: &Vec<T>) -> ExtendedMatrix<T, V>
    {
        let mut all_elements = Vec::new();
        for row_number in rows_numbers
        {
            let converted_row_number = conversion_uint_into_usize(*row_number);

            let node_dof_parameter =
                self.state.nodes_dof_parameters_global[converted_row_number];
            if let Some(position) = self.boundary_conditions
                .iter()
                .position(|bc|
                    bc.is_dof_parameter_data_same(
                        node_dof_parameter.dof_parameter(),
                        node_dof_parameter.node_number()))
            {
                let value = self.boundary_conditions[position].extract_value();
                all_elements.push(value);
            }
            else
            {
                all_elements.push(V::from(0f32));
            }
        }

        let mut converted_rows_numbers = T::from(0u8);
        (0..rows_numbers.len()).for_each(|_| converted_rows_numbers += T::from(1u8));

        let matrix = ExtendedMatrix::create(
            converted_rows_numbers,
            T::from(1u8),
            all_elements, self.state.tolerance);
        matrix
    }


    fn compose_displacements_matrix(&self, ua_matrix: ExtendedMatrix<T, V>,
        ub_matrix: ExtendedMatrix<T, V>, ua_ra_rows_numbers: &Vec<T>,
        ub_rb_rows_numbers: &Vec<T>) -> ExtendedMatrix<T, V>
    {
        let ua_values = ua_matrix.extract_all_elements_values();
        let ub_values = ub_matrix.extract_all_elements_values();
        let mut all_displacements_values =
            vec![V::from(0f32); self.state.nodes_dof_parameters_global.len()];

        let mut i = T::from(0u8);
        (0..ua_ra_rows_numbers.len()).for_each(|index|
            {
                let displacement_value = extract_element_value(
                    i, T::from(0u8), &ua_values);
                let converted_index = conversion_uint_into_usize(ua_ra_rows_numbers[index]);
                all_displacements_values[converted_index] = displacement_value;
                i += T::from(1u8);
            });

        let mut j = T::from(0u8);
        (0..ub_rb_rows_numbers.len()).for_each(|index|
            {
                let displacement_value = extract_element_value(
                    j, T::from(0u8), &ub_values);
                let converted_index = conversion_uint_into_usize(ub_rb_rows_numbers[index]);
                all_displacements_values[converted_index] = displacement_value;
                j += T::from(1u8);
            });

        let mut rows_number = T::from(0u8);
        (0..self.state.nodes_dof_parameters_global.len()).for_each(|_| rows_number += T::from(1u8));

        let displacement_matrix =
            ExtendedMatrix::create(
                rows_number, T::from(1u8), all_displacements_values,
                self.state.tolerance);
        displacement_matrix
    }


    pub fn global_analysis(&mut self) -> Result<GlobalAnalysisResult<T, V>, String>
    {
        self.update_nodes_dof_parameters_global()?;

        if self.boundary_conditions.iter().position(|bc|
            bc.is_type_same(BCType::Displacement)).is_none()
        {
            return Err("FEModel: Model could not be analyzed because there are no restraints were \
                applied!".into())
        }
        let mut global_stiffness_matrix =
            self.compose_global_stiffness_matrix()?;

        let removed_zeros_rows_columns =
            global_stiffness_matrix.remove_zeros_rows_columns();
        self.shrink_of_nodes_dof_parameters(&removed_zeros_rows_columns)?;
        let mut ub_rb_rows_numbers = Vec::new();
        let mut separation_positions = Vec::new();
        self.compose_separation_positions(&mut ub_rb_rows_numbers, &mut separation_positions);
        let mut ua_ra_rows_numbers = Vec::new();
        self.compose_ua_ra_rows_numbers(&ub_rb_rows_numbers, &mut ua_ra_rows_numbers);
        let ra_matrix = self.compose_matrix_by_rows_numbers(&ua_ra_rows_numbers);
        let ub_matrix = self.compose_matrix_by_rows_numbers(&ub_rb_rows_numbers);
        let separated_matrix =
            separate(global_stiffness_matrix, separation_positions, self.state.tolerance)?;
        let ua_matrix = separated_matrix.k_aa()
            .naive_gauss_elimination(&ra_matrix.add_subtract_matrix(
            &separated_matrix.k_ab().multiply_by_matrix(&ub_matrix)?,
            Operation::Subtraction)?)?;
        let reactions_values_matrix = separated_matrix.k_ba()
            .multiply_by_matrix(&ua_matrix)?
            .add_subtract_matrix(
                &separated_matrix.k_bb()
                    .multiply_by_matrix(&ub_matrix)?,
                        Operation::Addition)?;
        let all_reactions =
            reactions_values_matrix.extract_all_elements_values();
        let reactions_values_matrix_shape = reactions_values_matrix.get_shape();
        let mut reactions_values = Vec::new();

        let mut row = T::from(0u8);
        while row < reactions_values_matrix_shape.0
        {
            let mut column = T::from(0u8);
            while column < reactions_values_matrix_shape.1
            {
                let reaction_value = extract_element_value(row, column,
                    &all_reactions);
                reactions_values.push(reaction_value);
                column += T::from(1u8);
            }
            row += T::from(1u8);
        }

        let mut reactions_dof_parameters_data = Vec::new();
        for row_number in &ub_rb_rows_numbers
        {
            let converted_row_number = conversion_uint_into_usize(*row_number);
            reactions_dof_parameters_data.push(
                self.state.nodes_dof_parameters_global[converted_row_number]);
        }
        let displacements_dof_parameters_data =
            self.state.nodes_dof_parameters_global.clone();
        let displacements_values_matrix = self.compose_displacements_matrix(
            ua_matrix, ub_matrix, &ua_ra_rows_numbers, &ub_rb_rows_numbers);
        let all_displacements =
            displacements_values_matrix.extract_all_elements_values();
        let displacements_values_matrix_shape = displacements_values_matrix.get_shape();
        let mut displacements_values = Vec::new();

        let mut row = T::from(0u8);
        while row < displacements_values_matrix_shape.0
        {
            let mut column = T::from(0u8);
            while column < displacements_values_matrix_shape.1
            {
                let displacement_value = extract_element_value(row, column,
                    &all_displacements);
                displacements_values.push(displacement_value);
                column += T::from(1u8);
            }
            row += T::from(1u8);
        }

        let global_analysis_result =
            GlobalAnalysisResult::create(
                reactions_values, reactions_dof_parameters_data,
                displacements_values, displacements_dof_parameters_data);
        Ok(global_analysis_result)
    }


    pub fn elements_analysis(&self, global_displacements: &Displacements<T, V>)
        -> Result<ElementsAnalysisResult<T, V>, String>
    {
        let mut elements_analysis_result = HashMap::new();
        let mut analyzed_elements_types: HashMap<FEType, Vec<T>> = HashMap::new();
        for (element_number, element) in self.elements.iter()
        {
            let element_type = element.extract_fe_type();
            if let Some(element_numbers) = analyzed_elements_types
                .get_mut(&element_type)
            {
                element_numbers.push(*element_number);
            }
            else
            {
                analyzed_elements_types.insert(element_type, vec![*element_number]);
            }
            let element_analysis_data = element.extract_element_analysis_data(
                global_displacements, self.state.tolerance, &self.nodes)?;
            elements_analysis_result.insert(*element_number, element_analysis_data);
        }

        let elements_analysis_result = ElementsAnalysisResult::create(
            elements_analysis_result, analyzed_elements_types);

        Ok(elements_analysis_result)
    }


    pub fn beam_elements_nodal_forces(&self, beam_element_numbers: &[T],
        elements_analysis_data: &HashMap<T, ElementAnalysisData<V>>)
        -> HashMap<T, BeamElementNodalForces<T, V>>
    {
        let mut beam_elements_nodal_forces: HashMap<T, BeamElementNodalForces<T, V>> = HashMap::new();

        for beam_element_number in beam_element_numbers
        {
            let forces = elements_analysis_data
                .get(beam_element_number).unwrap()
                .forces_values().unwrap().to_vec();

            let moment_y_average = forces[4];
            let moment_y_min = forces[5];
            let moment_y_max = forces[6];
            let moment_z_average = forces[7];
            let moment_z_min = forces[8];
            let moment_z_max = forces[9];

            let nodes_numbers = self.elements
                .get(beam_element_number).unwrap()
                .extract_nodes_numbers();

            let node_1_number = nodes_numbers[0];
            let node_2_number = nodes_numbers[1];

            let node_1_forces_components =
                vec![ForceComponent::MomentY, ForceComponent::MomentZ];
            let node_2_forces_components =
                vec![ForceComponent::MomentY, ForceComponent::MomentZ];
            let mut node_1_forces_values = Vec::new();
            let mut node_2_forces_values = Vec::new();

            let mut adjacent_elements_numbers = Vec::new();
            let mut beam_elements_numbers_for_search = beam_element_numbers.to_vec();
            while let Some(position) = beam_elements_numbers_for_search.iter()
                .position(|number| number != beam_element_number &&
                    self.elements.get(number).unwrap()
                        .is_node_belong_element(node_1_number))
            {
                adjacent_elements_numbers.push(
                    beam_elements_numbers_for_search.remove(position));
            }

            if adjacent_elements_numbers.len() == 1
            {
                let adjacent_element_number = adjacent_elements_numbers[0];
                let adjacent_forces = elements_analysis_data
                    .get(&adjacent_element_number).unwrap()
                    .forces_values().unwrap().to_vec();
                let adjacent_moment_y_average = adjacent_forces[4];
                let adjacent_moment_z_average = adjacent_forces[7];

                if moment_y_average > adjacent_moment_y_average
                {
                    node_1_forces_values.push(moment_y_min);
                    node_2_forces_values.push(moment_y_max);
                }
                else
                {
                    node_1_forces_values.push(moment_y_max);
                    node_2_forces_values.push(moment_y_min);
                }

                if moment_z_average > adjacent_moment_z_average
                {
                    node_1_forces_values.push(moment_z_min);
                    node_2_forces_values.push(moment_z_max);
                }
                else
                {
                    node_1_forces_values.push(moment_z_max);
                    node_2_forces_values.push(moment_z_min);
                }

                let mut beam_element_nodal_forces = HashMap::new();

                beam_element_nodal_forces.insert(node_1_number, NodalForces::create(
                    node_1_forces_values, node_1_forces_components));
                beam_element_nodal_forces.insert(node_2_number, NodalForces::create(
                    node_2_forces_values, node_2_forces_components));

                beam_elements_nodal_forces.insert(*beam_element_number,
                    BeamElementNodalForces::create(beam_element_nodal_forces));

                continue;
            }

            let mut adjacent_elements_numbers = Vec::new();
            let mut beam_elements_numbers_for_search = beam_element_numbers.to_vec();
            while let Some(position) = beam_elements_numbers_for_search.iter()
                .position(|number| number != beam_element_number &&
                    self.elements.get(number).unwrap()
                        .is_node_belong_element(node_2_number))
            {
                adjacent_elements_numbers.push(
                    beam_elements_numbers_for_search.remove(position));
            }
            if adjacent_elements_numbers.len() == 1
            {
                let adjacent_element_number = adjacent_elements_numbers[0];
                let adjacent_forces = elements_analysis_data
                    .get(&adjacent_element_number).unwrap()
                    .forces_values().unwrap().to_vec();
                let adjacent_moment_y_average = adjacent_forces[4];
                let adjacent_moment_z_average = adjacent_forces[7];

                if moment_y_average > adjacent_moment_y_average
                {
                    node_2_forces_values.push(moment_y_min);
                    node_1_forces_values.push(moment_y_max);
                }
                else
                {
                    node_2_forces_values.push(moment_y_max);
                    node_1_forces_values.push(moment_y_min);
                }

                if moment_z_average > adjacent_moment_z_average
                {
                    node_2_forces_values.push(moment_z_min);
                    node_1_forces_values.push(moment_z_max);
                }
                else
                {
                    node_2_forces_values.push(moment_z_max);
                    node_1_forces_values.push(moment_z_min);
                }

                let mut beam_element_nodal_forces = HashMap::new();

                beam_element_nodal_forces.insert(node_1_number, NodalForces::create(
                    node_1_forces_values, node_1_forces_components));
                beam_element_nodal_forces.insert(node_2_number, NodalForces::create(
                    node_2_forces_values, node_2_forces_components));

                beam_elements_nodal_forces.insert(*beam_element_number,
                    BeamElementNodalForces::create(beam_element_nodal_forces));

                continue;
            }
        }

        beam_elements_nodal_forces
    }
}
