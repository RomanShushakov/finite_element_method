use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign, DivAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::collections::{HashSet, HashMap};
use std::iter::FromIterator;

use extended_matrix::extended_matrix::ExtendedMatrix;
use extended_matrix::matrix_element_position::MatrixElementPosition;
use extended_matrix::extended_matrix::Operation;
use extended_matrix::functions::{conversion_uint_into_usize, matrix_element_value_extractor};

use crate::fem::finite_elements::fe_node::{FENode, DeletedFENodeData};
use crate::fem::finite_elements::finite_element::{FiniteElement, FEType, DeletedFEData};
use crate::fem::global_analysis::fe_stiffness::{StiffnessGroup, StiffnessGroupKey};
use crate::fem::global_analysis::fe_boundary_condition::{BoundaryCondition, BCType, DeletedBCData};
use crate::fem::global_analysis::fe_global_analysis_result::{GlobalAnalysisResult, Displacements};
use crate::fem::global_analysis::fe_dof_parameter_data::
{
    global_dof, DOFParameterData, GLOBAL_DOF, GlobalDOFParameter
};

use crate::fem::element_analysis::fe_element_analysis_result::{ElementsAnalysisResult, EARType};

use crate::fem::functions::{separate, add_new_stiffness_sub_groups};

use crate::my_float::MyFloatTrait;


struct State<T, V>
{
    stiffness_groups: HashMap<StiffnessGroupKey<T>, Vec<MatrixElementPosition<T>>>,
    nodes_dof_parameters_global: Vec<DOFParameterData<T>>,
    tolerance: V,
}


impl<T, V> State<T, V>
{
    fn create(stiffness_groups: HashMap<StiffnessGroupKey<T>, Vec<MatrixElementPosition<T>>>,
              nodes_dof_parameters_global: Vec<DOFParameterData<T>>, tolerance: V,) -> Self
    {
        State { stiffness_groups, nodes_dof_parameters_global, tolerance }
    }
}


pub struct FEModel<T, V>
{
    nodes: HashMap<T, FENode<V>>,                           // Hashmap { node_number: Node }
    elements: HashMap<T, FiniteElement<T, V>>,              // Hashmap { element_number: FiniteElement }
    boundary_conditions: Vec<BoundaryCondition<T, V>>,
    state: State<T, V>,
}


impl<T, V> FEModel<T, V>
    where T: Copy + PartialEq + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + Eq + Hash +
             SubAssign + Debug + Mul<Output = T> + PartialOrd + Add<Output = T> + From<u8> +
             AddAssign + Ord + 'static,
          V: Copy + Sub<Output = V> + Mul<Output = V> + Add<Output = V> + Div<Output = V> +
             PartialEq + Debug + AddAssign + MulAssign + SubAssign + Into<f64> + PartialOrd +
             MyFloatTrait + From<f32> + MyFloatTrait<Other = V> + DivAssign + 'static,
{
    pub fn create(tolerance: V) -> Self
    {
        let state = State::create(HashMap::new(),
            Vec::new(), tolerance);

        FEModel { nodes: HashMap::new(), elements: HashMap::new(),
            boundary_conditions: Vec::new(), state }
    }


    pub fn reset(&mut self)
    {
        self.nodes = HashMap::new();
        self.elements = HashMap::new();
        self.boundary_conditions = Vec::new();
        self.state.stiffness_groups = HashMap::new();
        self.state.nodes_dof_parameters_global = Vec::new();
    }


    fn update_stiffness_groups(&mut self) -> Result<(), &str>
    {
        let mut stiffness_groups =
            HashMap::new();
        if self.nodes.len() < 2
        {
            self.state.stiffness_groups = HashMap::new();
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
                        add_new_stiffness_sub_groups(&mut stiffness_groups,
                            position, columns_number,
                            nodes_numbers[j], nodes_numbers[j])?;
                        position += T::from(1u8);
                    }
                    add_new_stiffness_sub_groups(&mut stiffness_groups,
                        position, columns_number,
                        excluded, v_lhs[j])?;
                    position += T::from(1u8);
                }
            }

            for i in 0..nodes_numbers.len() - 1
            {
                add_new_stiffness_sub_groups(&mut stiffness_groups,
                    position, columns_number,
                    nodes_numbers[nodes_numbers.len() - 1],
                        nodes_numbers[i])?;
                position += T::from(1u8);
            }
            add_new_stiffness_sub_groups(&mut stiffness_groups,
                position, columns_number,
                nodes_numbers[nodes_numbers.len() - 1],
                nodes_numbers[nodes_numbers.len() - 1])?;
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


    pub fn update_node(&mut self, number: T, x: V, y: V, z: V, optional_properties: Option<Vec<V>>)
        -> Result<(), String>
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
                .filter(|element|
                    element.is_node_belong_element(number))
            {
                if let Some(properties) = optional_properties.as_ref()
                {
                    match element.copy_fe_type()
                    {
                        FEType::Beam2n1ipT =>
                            {
                                if properties.len() != 11
                                {
                                    return Err("FEModel: Update node action: Incorrect length of \
                                        properties data!".into());
                                }
                                for (i, value) in properties.iter().enumerate()
                                {
                                    if *value <= V::from(0f32) && [i != 5, i < 8].iter()
                                        .all(|condition| *condition == true)
                                    {
                                        return Err("FEModel: Update node action: All properties \
                                            values should be greater than zero!".into());
                                    }
                                }
                            },
                        _ =>
                            {
                                if properties.len() < 2 || properties.len() > 3
                                {
                                    return Err("FEModel: Update node action: Incorrect \
                                        length of properties data!".into());
                                }

                                for value in properties.iter()
                                {
                                    if *value <= V::from(0f32)
                                    {
                                        return Err("FEModel: Update node action: All properties \
                                            values should be greater than zero!".into());
                                    }
                                }
                            },
                    }

                    let element_nodes_numbers = element.copy_nodes_numbers();
                    element.update(element_nodes_numbers, properties.clone(),
                        self.state.tolerance, &self.nodes)?;
                }
                else
                {
                    element.refresh(self.state.tolerance, &self.nodes)?;
                }
            }
            Ok(())
        }
        else
        {
            return Err(format!("FEModel: Node {:?} could not be updated because it does not \
                exist!", number));
        }
    }


    pub fn delete_node(&mut self, number: T) -> Result<(DeletedFENodeData<T, V>,
        Option<Vec<DeletedFEData<T, V>>>, Option<Vec<DeletedBCData<T, V>>>), String>
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

        let mut deleted_finite_elements_data = Vec::new();
        for element_number in elements_numbers_for_delete
        {
            let deleted_element = self.elements.remove(&element_number).unwrap();
            let deleted_element_number = element_number;
            let deleted_finite_element_data =
                DeletedFEData::create(deleted_element_number, deleted_element);
            deleted_finite_elements_data.push(deleted_finite_element_data);
        }

        let mut deleted_bcs_data = Vec::new();
        while let Some(position) = self.boundary_conditions.iter().position(|bc|
            bc.is_node_number_same(number))
        {
            let deleted_bc = self.boundary_conditions.remove(position);
            let deleted_bc_data = DeletedBCData::create(deleted_bc);
            deleted_bcs_data.push(deleted_bc_data);
        }

        let deleted_node = self.nodes.remove(&number).unwrap();
        self.update_stiffness_groups()?;

        let deleted_fe_node_data =
            DeletedFENodeData::create(number, deleted_node);

        let optional_deleted_finite_elements_data =
            {
                if deleted_finite_elements_data.is_empty()
                {
                    None
                }
                else
                {
                    Some(deleted_finite_elements_data)
                }
            };

        let optional_deleted_bcs_data =
            {
                if deleted_bcs_data.is_empty()
                {
                    None
                }
                else
                {
                    Some(deleted_bcs_data)
                }
            };

        Ok((deleted_fe_node_data, optional_deleted_finite_elements_data, optional_deleted_bcs_data))
    }


    pub fn add_element(&mut self, element_number: T, element_type: FEType, nodes_numbers: Vec<T>,
        properties: Vec<V>) -> Result<(), String>
    {
        if self.elements.contains_key(&element_number)
        {
            return Err(format!("FEModel: Element {:?} could not be added! The element with the \
                same number does already exist!", element_number));
        }

        let nodes_numbers_set = HashSet::<T>::from_iter(
            nodes_numbers.iter().cloned());
        if nodes_numbers.len() != nodes_numbers_set.len()
        {
            return Err(format!("FEModel: Element {:?} could not be added! All nodes numbers \
                should be unique!", element_number));
        }

        match element_type
        {
            FEType::Beam2n1ipT =>
                {
                    if nodes_numbers.len() != 2
                    {
                        return Err(format!("FEModel: Element {:?} could not be added! \
                            Incorrect number of nodes!", element_number));
                    }

                    if properties.len() != 11
                    {
                        return Err(format!("FEModel: Element {:?} could not be added! \
                            Incorrect length of properties data!", element_number));
                    }
                    for (i, value) in properties.iter().enumerate()
                    {
                        if *value <= V::from(0f32) && [i != 5, i < 8].iter()
                            .all(|condition| *condition == true)
                        {
                            return Err(format!("FEModel: All properties values for element {:?} \
                                should be greater than zero!", element_number));
                        }
                    }
                },
            _ =>
                {
                    if nodes_numbers.len() != 2
                    {
                        return Err(format!("FEModel: Element {:?} could not be added! \
                            Incorrect number of nodes!", element_number));
                    }

                    if properties.len() < 2 || properties.len() > 3
                    {
                        return Err(format!("FEModel: Element {:?} could not be added! \
                            Incorrect length of properties data!", element_number));
                    }

                    for value in properties.iter()
                    {
                        if *value <= V::from(0f32)
                        {
                            return Err(format!("FEModel: All properties values for element {:?} \
                                should be greater than zero!", element_number));
                        }
                    }
                }
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

        if self.elements.iter().position(|(number, element)|
            element.is_nodes_numbers_same(nodes_numbers.clone()) &&
                *number != element_number).is_some()
        {
            return Err(format!("FEModel: Element {:?} could not be updated! The element with \
                    the same nodes numbers does already exist!", element_number));
        }

        if let Some(element) = self.elements.get_mut(&element_number)
        {
            match element.copy_fe_type()
            {
                FEType::Beam2n1ipT =>
                    {
                        if nodes_numbers.len() != 2
                        {
                            return Err(format!("FEModel: Element {:?} could not be updated! \
                                Incorrect number of nodes!", element_number));
                        }

                        if properties.len() != 11
                        {
                            return Err(format!("FEModel: Element {:?} could not be updated! \
                                Incorrect length of properties data!", element_number));
                        }

                        for (i, value) in properties.iter().enumerate()
                        {
                            if *value <= V::from(0f32) && [i != 5, i < 8]
                                .iter()
                                .all(|condition| *condition == true)
                            {
                                return Err(format!("FEModel: All properties values for element {:?} \
                                    should be greater than zero!", element_number));
                            }
                        }
                    },
                _ =>
                    {
                        if nodes_numbers.len() != 2
                        {
                            return Err(format!("FEModel: Element {:?} could not be updated! \
                                Incorrect number of nodes!", element_number));
                        }

                        if properties.len() < 2 || properties.len() > 3
                        {
                             return Err(format!("FEModel: Element {:?} could not be updated! \
                                Incorrect length of properties data!", element_number));
                        }

                        for value in properties.iter()
                        {
                            if *value <= V::from(0f32)
                            {
                                return Err(format!("FEModel: All properties values for element {:?} \
                                    should be greater than zero!", element_number));
                            }
                        }
                    }
            }
            element.update(nodes_numbers, properties, self.state.tolerance, &self.nodes)?;
            Ok(())
        }
        else
        {
            return Err(format!("FEModel: Element {:?} could not be updated because it does not \
                exist!", element_number));
        }
    }


    pub fn delete_element(&mut self, number: T) -> Result<DeletedFEData<T, V>, String>
    {

        if let Some(deleted_element) = self.elements.remove(&number)
        {
            let deleted_element_number = number;
            Ok(DeletedFEData::create(deleted_element_number, deleted_element))
        }
        else
        {
            Err(format!("FEModel: Element {:?} could not be deleted because it does not \
                exist!", number))
        }
    }


    fn compose_global_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, String>
    {
        if self.elements.is_empty()
        {
            return Err("FEModel: Global stiffness matrix could not be composed because there are \
                no elements in the model!".to_string());
        }

        if self.nodes.keys().any(|node_number|
            self.elements.values().position(|element|
                element.is_node_belong_element(*node_number)).is_none())
        {
            return Err("FEModel: Global stiffness matrix could not be composed because there are \
                free nodes exist!".to_string());
        }

        let mut nodes_len_value = T::from(0u8);
        (0..self.nodes.len()).for_each(|_| nodes_len_value += T::from(1u8));

        let mut global_stiffness_matrix = ExtendedMatrix::<T, V>::create(
            nodes_len_value * global_dof::<T>(),
            nodes_len_value * global_dof::<T>(),
            vec![V::from(0f32); (self.nodes.len() * GLOBAL_DOF).pow(2)],
            self.state.tolerance)?;

        for element in self.elements.values()
        {
            let mut element_stiffness_matrix = element.extract_stiffness_matrix()?;

            let element_stiffness_groups = element.extract_stiffness_groups();

            for element_stiffness_group in element_stiffness_groups
            {
                let stiffness_group_key = StiffnessGroupKey {
                    stiffness_type: element_stiffness_group.stiffness_type,
                    number_1: element_stiffness_group.number_1,
                    number_2: element_stiffness_group.number_2 };
                if let Some(matrix_elements_positions) =
                    self.state.stiffness_groups.get(&stiffness_group_key)
                {
                    global_stiffness_matrix.add_submatrix_to_assemblage(
                        &mut element_stiffness_matrix,
                        matrix_elements_positions,
                        &element_stiffness_group.positions);
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


    pub fn delete_bc(&mut self, bc_type: BCType, number: T) -> Result<DeletedBCData<T, V>, String>
    {
        if let Some(position) =  self.boundary_conditions.iter().position(|bc|
            bc.is_number_same(number) && bc.is_type_same(bc_type))
        {
            let deleted_bc = self.boundary_conditions.remove(position);
            Ok(DeletedBCData::create(deleted_bc))
        }
        else
        {
            Err(format!("FEModel: {} could not be deleted because current {} number does not \
                exist!", bc_type.as_str(), bc_type.as_str().to_lowercase()))
        }
    }


    fn shrink_of_nodes_dof_parameters(&mut self, zeros_rows_columns: &Vec<MatrixElementPosition<T>>)
        -> Result<(), String>
    {
        for row_column in zeros_rows_columns
        {
            let mut row_column_as_index = 0usize;
            let mut n = T::from(0u8);
            while n < *row_column.ref_column()
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
                        dof_parameter_data.copy_dof_parameter(),
                        dof_parameter_data.copy_node_number()))
            {
                let bc_type = self.boundary_conditions[position].copy_bc_type();
                let dof_parameter = dof_parameter_data.copy_dof_parameter();
                let node_number = dof_parameter_data.copy_node_number();
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
                        dof_parameter_data.copy_dof_parameter(),
                        dof_parameter_data.copy_node_number())
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


    fn compose_ua_ra_rows_numbers(&self, ub_rb_rows_numbers: &[T],
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


    fn compose_matrix_by_rows_numbers(&self, rows_numbers: &[T], bc_type: BCType)
        -> Result<ExtendedMatrix<T, V>, String>
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
                        node_dof_parameter.copy_dof_parameter(),
                        node_dof_parameter.copy_node_number()) &&
                        bc.is_type_same(bc_type)
                )
            {
                let value = self.boundary_conditions[position].copy_value();
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
            all_elements, self.state.tolerance)?;
        Ok(matrix)
    }


    fn compose_displacements_matrix(&self, ua_matrix: ExtendedMatrix<T, V>,
        ub_matrix: ExtendedMatrix<T, V>, ua_ra_rows_numbers: &[T],
        ub_rb_rows_numbers: &Vec<T>) -> Result<ExtendedMatrix<T, V>, String>
    {
        let mut all_displacements_values =
            vec![V::from(0f32); self.state.nodes_dof_parameters_global.len()];

        let mut i = T::from(0u8);
        let mut index = 0usize;
        while index < ua_ra_rows_numbers.len()
        {
            let displacement_value = matrix_element_value_extractor(i, T::from(0u8), &ua_matrix)?;
            let converted_index = conversion_uint_into_usize(ua_ra_rows_numbers[index]);
            all_displacements_values[converted_index] = displacement_value;
            i += T::from(1u8);
            index += 1usize;
        }

        let mut j = T::from(0u8);
        let mut index = 0usize;
        while index < ub_rb_rows_numbers.len()
        {
            let displacement_value = matrix_element_value_extractor(j, T::from(0u8), &ub_matrix)?;
            let converted_index = conversion_uint_into_usize(ub_rb_rows_numbers[index]);
            all_displacements_values[converted_index] = displacement_value;
            j += T::from(1u8);
            index += 1usize;
        }

        let mut rows_number = T::from(0u8);
        (0..self.state.nodes_dof_parameters_global.len()).for_each(|_| rows_number += T::from(1u8));

        let displacement_matrix =
            ExtendedMatrix::create(
                rows_number, T::from(1u8), all_displacements_values,
                self.state.tolerance);
        displacement_matrix
    }


    pub fn global_analysis(&mut self , colsol_usage: bool)
        -> Result<GlobalAnalysisResult<T, V>, String>
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

        let ra_matrix = self.compose_matrix_by_rows_numbers(
            &ua_ra_rows_numbers, BCType::Force)?;
        let ub_matrix = self.compose_matrix_by_rows_numbers(
            &ub_rb_rows_numbers, BCType::Displacement)?;
        let rb_c_matrix = self.compose_matrix_by_rows_numbers(
            &ub_rb_rows_numbers, BCType::Force)?;

        let separated_matrix =
            separate(global_stiffness_matrix, separation_positions, self.state.tolerance)?;

        let ua_matrix = separated_matrix.ref_k_aa()
            .direct_solution(&ra_matrix.add_subtract_matrix(
                &separated_matrix.ref_k_ab().multiply_by_matrix(&ub_matrix)?,
                Operation::Subtraction)?, colsol_usage)?;

        let reactions_values_matrix = separated_matrix.ref_k_ba()
            .multiply_by_matrix(&ua_matrix)?
            .add_subtract_matrix(
                &separated_matrix.ref_k_bb()
                    .multiply_by_matrix(&ub_matrix)?,
                        Operation::Addition)?
            .add_subtract_matrix(&rb_c_matrix, Operation::Subtraction)?;

        let reactions_values_matrix_shape = reactions_values_matrix.copy_shape();
        let mut reactions_values = Vec::new();

        let mut row = T::from(0u8);
        while row < reactions_values_matrix_shape.0
        {
            let mut column = T::from(0u8);
            while column < reactions_values_matrix_shape.1
            {
                let reaction_value = matrix_element_value_extractor(row, column, &reactions_values_matrix)?;
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
            ua_matrix, ub_matrix, &ua_ra_rows_numbers, &ub_rb_rows_numbers)?;

        let displacements_values_matrix_shape = displacements_values_matrix.copy_shape();
        let mut displacements_values = Vec::new();

        let mut row = T::from(0u8);
        while row < displacements_values_matrix_shape.0
        {
            let mut column = T::from(0u8);
            while column < displacements_values_matrix_shape.1
            {
                let displacement_value = matrix_element_value_extractor(row, column, &displacements_values_matrix)?;
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
        let mut elements_analysis_result = ElementsAnalysisResult::create();
        for (element_number, element) in self.elements.iter()
        {
            let fe_type = element.copy_fe_type();
            let element_analysis_data = element.extract_element_analysis_data(
                global_displacements, self.state.tolerance, &self.nodes)?;

            elements_analysis_result.add_to_analysis_data(*element_number, element_analysis_data);
            elements_analysis_result.add_to_types(fe_type, *element_number)?;
        }

        Ok(elements_analysis_result)
    }


    pub fn copy_node_coordinates(&self, node_number: &T) -> Result<(V, V, V), String>
    {
        if let Some(node) = self.nodes.get(node_number)
        {
            Ok(node.copy_coordinates())
        }
        else
        {
            Err(format!("FEModel: Node with number {:?} does not exist!", node_number))
        }
    }


    pub fn copy_element_type(&self, element_number: &T) -> Result<FEType, String>
    {
        if let Some(element) = self.elements.get(element_number)
        {
            Ok(element.copy_fe_type())
        }
        else
        {
            Err(format!("FEModel: Element with number {:?} does not exist!", element_number))
        }
    }


    pub fn copy_element_nodes_numbers(&self, element_number: &T) -> Result<Vec<T>, String>
    {
        if let Some(element) = self.elements.get(element_number)
        {
            Ok(element.copy_nodes_numbers())
        }
        else
        {
            Err(format!("FEModel: Element with number {:?} does not exist!", element_number))
        }
    }


    pub fn copy_element_properties(&self, element_number: &T) -> Result<Vec<V>, String>
    {
        if let Some(element) = self.elements.get(element_number)
        {
            Ok(element.copy_properties())
        }
        else
        {
            Err(format!("FEModel: Element with number {:?} does not exist!", element_number))
        }
    }


    pub fn extract_unique_elements_of_rotation_matrix(&self, element_number: &T)
        -> Result<Vec<V>, String>
    {
        if let Some(element) = self.elements.get(element_number)
        {
            Ok(element.extract_unique_elements_of_rotation_matrix()?)
        }
        else
        {
            Err(format!("FEModel: Element with number {:?} does not exist!", element_number))
        }
    }


    pub fn copy_bc_node_number(&self, bc_type: BCType, number: T) -> Result<T, String>
    {
        if let Some(position) = self.boundary_conditions.iter()
            .position(|bc| bc.is_type_same(bc_type) && bc.
                is_number_same(number))
        {
            Ok(self.boundary_conditions[position].copy_node_number())
        }
        else
        {
            Err(format!("FEModel: {:?} boundary condition with number {:?} does not exist!",
                bc_type.as_str(), number))
        }
    }


    pub fn copy_bc_dof_parameter(&self, bc_type: BCType, number: T)
        -> Result<GlobalDOFParameter, String>
    {
        if let Some(position) = self.boundary_conditions.iter()
            .position(|bc| bc.is_type_same(bc_type) &&
                bc.is_number_same(number))
        {
            Ok(self.boundary_conditions[position].copy_dof_parameter())
        }
        else
        {
            Err(format!("FEModel: {:?} boundary condition with number {:?} does not exist!",
                bc_type.as_str(), number))
        }
    }


    pub fn copy_bc_value(&self, bc_type: BCType, number: T) -> Result<V, String>
    {
        if let Some(position) = self.boundary_conditions.iter()
            .position(|bc| bc.is_type_same(bc_type) && bc.is_number_same(number))
        {
            Ok(self.boundary_conditions[position].copy_value())
        }
        else
        {
            Err(format!("FEModel: {:?} boundary condition with number {:?} does not exist!",
                bc_type.as_str(), number))
        }
    }


    pub fn is_node_number_exist(&self, number: T) -> bool
    {
        self.nodes.contains_key(&number)
    }


    pub fn is_element_number_exist(&self, number: T) -> bool
    {
        self.elements.contains_key(&number)
    }


    pub fn is_bc_key_exist(&self, number: T, bc_type: BCType) -> bool
    {
        if self.boundary_conditions.iter().position(|bc|
            bc.is_number_same(number) && bc.is_type_same(bc_type)).is_some()
        {
            return true;
        }
        false
    }


    pub fn extract_all_nodes_numbers(&self) -> Vec<T>
    {
        let mut nodes_numbers = Vec::new();
        for node_number in self.nodes.keys()
        {
            nodes_numbers.push(*node_number);
        }
        nodes_numbers
    }


    pub fn extract_all_elements_numbers(&self) -> Vec<T>
    {
        let mut elements_numbers = Vec::new();
        for element_number in self.elements.keys()
        {
            elements_numbers.push(*element_number);
        }
        elements_numbers
    }


    pub fn extract_all_bc_types_numbers(&self) -> Vec<(BCType, T)>
    {
        let mut bc_types_numbers = Vec::new();
        for bc in self.boundary_conditions.iter()
        {
            let bc_type = bc.copy_bc_type();
            let bc_number = bc.copy_number();
            bc_types_numbers.push((bc_type, bc_number));
        }
        bc_types_numbers
    }
}
