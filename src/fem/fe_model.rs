use std::collections::{HashSet, HashMap};
use std::iter::FromIterator;

use extended_matrix::{Matrix, Position, FloatTrait, BasicOperationsTrait, TryIntoSymmetricCompactedMatrixTrait, Vector};
use colsol::{factorization, find_unknown};

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

use crate::fem::functions::
{
    separate, add_new_stiffness_sub_groups, is_points_of_quadrilateral_on_the_same_line,
    is_points_of_quadrilateral_on_the_same_plane, convex_hull_on_four_points_on_plane,
};

use crate::fem::separated_matrix::SeparatedMatrix;

use super::finite_elements::beam::beam2n1ipt::Beam2n1ipT;
use super::finite_elements::plate::plate4n4ip::Plate4n4ip;


struct State<V>
{
    stiffness_groups: HashMap<StiffnessGroupKey, Vec<Position>>,
    nodes_dof_parameters_global: Vec<DOFParameterData>,
    optional_ua_ra_rows_numbers: Option<Vec<u32>>,
    optional_ub_rb_rows_numbers: Option<Vec<u32>>,
    optional_ua_matrix: Option<Vector<V>>,
    optional_ub_matrix: Option<Vector<V>>,
    optional_ra_matrix: Option<Vector<V>>,
    optional_rb_c_matrix: Option<Vector<V>>,
    optional_separated_matrix: Option<SeparatedMatrix<V>>,
    rel_tol: V,
    abs_tol: V,
}


impl<V> State<V>
{
    fn create(
        stiffness_groups: HashMap<StiffnessGroupKey, Vec<Position>>,
        nodes_dof_parameters_global: Vec<DOFParameterData>, 
        rel_tol: V,
        abs_tol: V,
    ) 
        -> Self
    {
        let optional_ua_ra_rows_numbers = None;
        let optional_ub_rb_rows_numbers = None;
        let optional_ua_matrix = None;
        let optional_ub_matrix = None;
        let optional_ra_matrix = None;
        let optional_rb_c_matrix = None;
        let optional_separated_matrix = None;
        State 
        { 
            stiffness_groups, 
            nodes_dof_parameters_global, 
            optional_ua_ra_rows_numbers, 
            optional_ub_rb_rows_numbers, 
            optional_ua_matrix, 
            optional_ub_matrix, 
            optional_ra_matrix, 
            optional_rb_c_matrix, 
            optional_separated_matrix,
            rel_tol,
            abs_tol,
        }
    }
}


pub struct FEModel<V>
{
    nodes: HashMap<u32, FENode<V>>,                           // Hashmap { node_number: Node }
    elements: HashMap<u32, FiniteElement<V>>,              // Hashmap { element_number: FiniteElement }
    boundary_conditions: Vec<BoundaryCondition<V>>,
    state: State<V>,
}


impl<V> FEModel<V>
    where V: FloatTrait<Output = V>
{
    pub fn create(rel_tol: V, abs_tol: V) -> Self
    {
        let state = State::create(
            HashMap::new(), Vec::new(), rel_tol, abs_tol,
        );

        FEModel { nodes: HashMap::new(), elements: HashMap::new(), boundary_conditions: Vec::new(), state }
    }


    fn reset_optional_state_values(&mut self)
    {
        self.state.optional_ua_ra_rows_numbers = None;
        self.state.optional_ub_rb_rows_numbers = None;
        self.state.optional_ua_matrix = None;
        self.state.optional_ub_matrix = None;
        self.state.optional_ra_matrix = None;
        self.state.optional_rb_c_matrix = None;
        self.state.optional_separated_matrix = None;
    }


    pub fn reset(&mut self)
    {
        self.nodes = HashMap::new();
        self.elements = HashMap::new();
        self.boundary_conditions = Vec::new();
        self.state.stiffness_groups = HashMap::new();
        self.state.nodes_dof_parameters_global = Vec::new();
        self.reset_optional_state_values();
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
            let mut position = 0;

            let mut columns_number = nodes_numbers.len();

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
                        add_new_stiffness_sub_groups(
                            &mut stiffness_groups,
                            position, 
                            columns_number,
                            nodes_numbers[j], 
                            nodes_numbers[j],
                        )?;
                        position += 1;
                    }
                    add_new_stiffness_sub_groups(
                        &mut stiffness_groups,
                        position, 
                        columns_number,
                        excluded, 
                        v_lhs[j],
                    )?;
                    position += 1;
                }
            }

            for i in 0..nodes_numbers.len() - 1
            {
                add_new_stiffness_sub_groups(
                    &mut stiffness_groups,
                    position, 
                    columns_number,
                    nodes_numbers[nodes_numbers.len() - 1],
                    nodes_numbers[i],
                )?;
                position += 1;
            }
            add_new_stiffness_sub_groups(
                &mut stiffness_groups,
                position, 
                columns_number,
                nodes_numbers[nodes_numbers.len() - 1],
                nodes_numbers[nodes_numbers.len() - 1],
            )?;
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
                let dof_parameter = GlobalDOFParameter::iterator()
                    .nth(dof)
                    .ok_or("FEModel: Could not find dof parameter!")?;
                let dof_parameter_data = DOFParameterData::create(*node_number, *dof_parameter);
                nodes_dof_parameters.push(dof_parameter_data);
            }
        }
        self.state.nodes_dof_parameters_global = nodes_dof_parameters;
        Ok(())
    }


    pub fn add_node(&mut self, number: u32, x: V, y: V, z: V, stiffness_groups_update: bool) -> Result<(), String>
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

        if stiffness_groups_update
        {
            self.update_stiffness_groups()?;
        }
        Ok(())

    }


    pub fn update_node(&mut self, number: u32, x: V, y: V, z: V) -> Result<(), String>
    {
        if self.nodes
            .iter()
            .position(|(node_number, node)| 
                *node_number != number && 
                node.is_coordinates_same(x, y, z))
            .is_some()
        {
            return Err(format!("FEModel: Node {:?} could not be updated because the node with the \
                same coordinates does already exist!", number))
        }

        if let Some(node) = self.nodes.get_mut(&number)
        {
            node.update(x, y, z);
            Ok(())
        }
        else
        {
            return Err(format!("FEModel: Node {:?} could not be updated because it does not \
                exist!", number));
        }
    }


    pub fn delete_node(&mut self, number: u32, stiffness_groups_update: bool) 
        -> Result<(DeletedFENodeData<V>, Option<Vec<DeletedFEData<V>>>, Option<Vec<DeletedBCData<V>>>), String>
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
            let deleted_finite_element_data = DeletedFEData::create(
                deleted_element_number, deleted_element,
            );
            deleted_finite_elements_data.push(deleted_finite_element_data);
        }

        let mut deleted_bcs_data = Vec::new();
        while let Some(position) = self.boundary_conditions
            .iter()
            .position(|bc| bc.is_node_number_same(number))
        {
            let deleted_bc = self.boundary_conditions.remove(position);
            let deleted_bc_data = DeletedBCData::create(deleted_bc);
            deleted_bcs_data.push(deleted_bc_data);
        }

        let deleted_node = self.nodes.remove(&number).unwrap();

        if stiffness_groups_update
        {
            self.update_stiffness_groups()?;
        }

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


    pub fn add_element(
        &mut self, element_number: u32, element_type: FEType, nodes_numbers: Vec<u32>, properties: Vec<V>,
    ) 
        -> Result<(), String>
    {
        if self.elements.contains_key(&element_number)
        {
            return Err(format!("FEModel: Element {:?} could not be added! The element with the \
                same number does already exist!", element_number));
        }

        let nodes_numbers_set = HashSet::<u32>::from_iter(nodes_numbers.iter().cloned());
        if nodes_numbers.len() != nodes_numbers_set.len()
        {
            return Err(format!("FEModel: Element {:?} could not be added! All nodes numbers \
                should be unique!", element_number));
        }

        match element_type
        {
            FEType::Truss2n1ip | FEType::Truss2n2ip =>
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
                },
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
            FEType::Mem4n4ip => 
                {
                    if nodes_numbers.len() != 4
                    {
                        return Err(format!("FEModel: Element {:?} could not be added! \
                            Incorrect number of nodes!", element_number));
                    }

                    if properties.len() != 3
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
                },
            FEType::Plate4n4ip => 
                {
                    if nodes_numbers.len() != 4
                    {
                        return Err(format!("FEModel: Element {:?} could not be added! \
                            Incorrect number of nodes!", element_number));
                    }

                    if properties.len() != 4
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
                },
        }

        if self.elements
            .values()
            .position(|element| 
                element.is_type_same(&element_type) && 
                element.is_nodes_numbers_same(nodes_numbers.clone()))
            .is_some()
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

        match element_type
        {
            FEType::Mem4n4ip | FEType::Plate4n4ip => 
            {
                let node_1 = self.nodes.get(&nodes_numbers[0]).unwrap();
                let node_2 = self.nodes.get(&nodes_numbers[1]).unwrap();
                let node_3 = self.nodes.get(&nodes_numbers[2]).unwrap();
                let node_4 = self.nodes.get(&nodes_numbers[3]).unwrap();
                
                if is_points_of_quadrilateral_on_the_same_line(
                    &[node_1.copy_x(), node_1.copy_y(), node_1.copy_z()], 
                    &[node_2.copy_x(), node_2.copy_y(), node_2.copy_z()],  
                    &[node_3.copy_x(), node_3.copy_y(), node_3.copy_z()], 
                    &[node_4.copy_x(), node_4.copy_y(), node_4.copy_z()], 
                    self.state.tolerance)
                {
                    return Err(format!("FEModel: Element {:?} could not be added! Three of \
                        four nodes lie on the same line!", element_number));
                }

                if !is_points_of_quadrilateral_on_the_same_plane(
                    &[node_1.copy_x(), node_1.copy_y(), node_1.copy_z()], 
                    &[node_2.copy_x(), node_2.copy_y(), node_2.copy_z()],  
                    &[node_3.copy_x(), node_3.copy_y(), node_3.copy_z()], 
                    &[node_4.copy_x(), node_4.copy_y(), node_4.copy_z()], 
                    self.state.tolerance)
                {
                    return Err(format!("FEModel: Element {:?} could not be added! Not all nodes \
                        lie on the same plane!", element_number));
                }

                if convex_hull_on_four_points_on_plane(
                    &nodes_numbers, 
                    &[
                        &[node_1.copy_x(), node_1.copy_y(), node_1.copy_z()], 
                        &[node_2.copy_x(), node_2.copy_y(), node_2.copy_z()],  
                        &[node_3.copy_x(), node_3.copy_y(), node_3.copy_z()], 
                        &[node_4.copy_x(), node_4.copy_y(), node_4.copy_z()]
                    ], 
                    self.state.rel_tol,
                    self.state.abs_tol,
                )?
                .len() != 4
                {
                    return Err(format!("FEModel: Element {:?} could not be added! \
                        Quadrilateral is not convex!", element_number));
                }
            }, 
            _ => (),
        }

        let element = FiniteElement::create(
            element_type, nodes_numbers, properties, &self.nodes, self.state.rel_tol, self.state.abs_tol,
        )?;
        self.elements.insert(element_number, element);

        Ok(())
    }


    pub fn update_element(
        &mut self, element_number: u32, nodes_numbers: Vec<u32>, properties: Vec<V>,
    ) 
        -> Result<(), String>
    {
        let nodes_numbers_set = HashSet::<u32>::from_iter(nodes_numbers.iter().cloned());
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

        if self.elements
            .iter()
            .position(|(number, element)| 
                element.is_nodes_numbers_same(nodes_numbers.clone()) && *number != element_number)
            .is_some()
        {
            return Err(format!("FEModel: Element {:?} could not be updated! The element with \
                the same nodes numbers does already exist!", element_number));
        }

        if let Some(element) = self.elements.get_mut(&element_number)
        {
            match element.copy_fe_type()
            {
                FEType::Truss2n1ip | FEType::Truss2n2ip =>
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
                },
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
                FEType::Mem4n4ip => 
                {
                    if nodes_numbers.len() != 4
                    {
                        return Err(format!("FEModel: Element {:?} could not be updated! \
                            Incorrect number of nodes!", element_number));
                    }

                    if properties.len() != 3
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
                },
            FEType::Plate4n4ip => 
                {
                    if nodes_numbers.len() != 4
                    {
                        return Err(format!("FEModel: Element {:?} could not be updated! \
                            Incorrect number of nodes!", element_number));
                    }

                    if properties.len() != 4
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
                },

            }
            element.update(nodes_numbers, properties, &self.nodes, self.state.rel_tol, self.state.abs_tol)?;
            Ok(())
        }
        else
        {
            return Err(format!("FEModel: Element {:?} could not be updated because it does not \
                exist!", element_number));
        }
    }


    pub fn delete_element(&mut self, number: u32) -> Result<DeletedFEData<V>, String>
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


    fn compose_global_stiffness_matrix(&self) -> Result<(Matrix<V>, Vec<u32>, Vec<u32>), String>
    {
        if self.elements.is_empty()
        {
            return Err("FEModel: Global stiffness matrix could not be composed because there are \
                no elements in the model!".to_string());
        }

        if self.nodes.keys().any(|node_number| self.elements
            .values()
            .position(|element| element.is_node_belong_element(*node_number))
            .is_none())
        {
            return Err("FEModel: Global stiffness matrix could not be composed because there are \
                free nodes exist!".to_string());
        }

        let mut nodes_len_value = self.nodes.len();

        let mut row_column_number = 0;
        let mut zero_rows_numbers = Vec::new();
        let mut zero_columns_numbers = Vec::new();
        while row_column_number < nodes_len_value * global_dof()
        {
            zero_rows_numbers.push(row_column_number);
            zero_columns_numbers.push(row_column_number);
            row_column_number += 1;
        }

        let mut global_stiffness_matrix = Matrix::<V>::create(
            nodes_len_value * global_dof(),
            nodes_len_value * global_dof(),
            &vec![V::from(0f32); (self.nodes.len() * GLOBAL_DOF).pow(2)],
        );

        for element in self.elements.values()
        {
            let mut element_stiffness_matrix = element.extract_stiffness_matrix()?;

            let element_stiffness_groups = element.extract_stiffness_groups();

            for element_stiffness_group in element_stiffness_groups
            {
                let stiffness_group_key = StiffnessGroupKey {
                    stiffness_type: element_stiffness_group.stiffness_type,
                    number_1: element_stiffness_group.number_1,
                    number_2: element_stiffness_group.number_2 
                };
                let element_matrix_elements_positions = element_stiffness_group.positions.clone();
                if let Some(matrix_elements_positions) = self.state.stiffness_groups
                    .get(&stiffness_group_key)
                {
                    for (matrix_element_position, element_matrix_element_position) in 
                        matrix_elements_positions.iter().zip(element_matrix_elements_positions.into_iter())
                    {
                        let element_value = element_stiffness_matrix.get_element_value(
                            &element_matrix_element_position,
                        )?;

                        let nonzero_row_number = matrix_element_position.0;
                        if let Some(position) = zero_rows_numbers
                            .iter()
                            .position(|row_number| 
                                row_number == nonzero_row_number && element_value != V::from(0f32))
                        {
                            zero_rows_numbers.remove(position);
                        }

                        let nonzero_column_number = matrix_element_position.1;
                        if let Some(position) = zero_columns_numbers
                            .iter()
                            .position(|column_number| 
                                column_number == nonzero_column_number && element_value != V::from(0f32))
                        {
                            zero_columns_numbers.remove(position);
                        }
                    }

                    for (lhs_position, rhs_position) in
                        matrix_elements_positions.iter().zip(element_stiffness_group.positions.iter())
                    {
                        let mut lhs_element_value = global_stiffness_matrix.get_mut_element_value(lhs_position)?;
                        let rhs_element_value = element_stiffness_matrix.get_element_value(rhs_position)?;
                        *lhs_element_value += rhs_element_value;
                    }
                }
            }
        }

        Ok((global_stiffness_matrix, zero_rows_numbers, zero_columns_numbers))
    }


    pub fn add_bc(
        &mut self, 
        bc_type: BCType, 
        number: u32, 
        node_number: u32,
        dof_parameter: GlobalDOFParameter, 
        value: V,
    ) 
        -> Result<(), String>
    {
        if self.boundary_conditions
            .iter()
            .position(|bc| bc.is_number_same(number) && bc.is_type_same(bc_type))
            .is_some()
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


    pub fn update_bc(
        &mut self, 
        bc_type: BCType, 
        number: u32, 
        node_number: u32,
        dof_parameter: GlobalDOFParameter, 
        value: V,
    ) 
        -> Result<(), String>
    {
        if !self.nodes.contains_key(&node_number)
        {
            return Err(format!("FEModel: {} could not be updated because the current node number \
                does not exist!", bc_type.as_str()));
        }

        if let Some(position) =  self.boundary_conditions
            .iter()
            .position(|bc| bc.is_number_same(number) && bc.is_type_same(bc_type))
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


    pub fn delete_bc(&mut self, bc_type: BCType, number: u32) -> Result<DeletedBCData<V>, String>
    {
        if let Some(position) =  self.boundary_conditions
            .iter()
            .position(|bc| bc.is_number_same(number) && bc.is_type_same(bc_type))
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


    fn shrink_of_nodes_dof_parameters(&mut self, zeros_rows_columns: &Vec<Position>) -> Result<(), String>
    {
        for row_column in zeros_rows_columns
        {
            let mut row_column_as_index = *row_column.1;

            let dof_parameter_data = self.state.nodes_dof_parameters_global
                .remove(row_column_as_index);

            if let Some(position) = self.boundary_conditions
                .iter()
                .position(|bc| bc
                    .is_dof_parameter_data_same(
                        dof_parameter_data.copy_dof_parameter(),
                        dof_parameter_data.copy_node_number()),
                    )
            {
                let bc_type = self.boundary_conditions[position].copy_bc_type();
                let dof_parameter = dof_parameter_data.copy_dof_parameter();
                let node_number = dof_parameter_data.copy_node_number();
                return Err(format!("FEModel: Model could not be analyzed because where are \
                    no stiffness to withstand {}::{:?} applied at node {:?}!", 
                        bc_type.as_str(),
                        dof_parameter, 
                        node_number))
            }
        }
        Ok(())
    }


    fn compose_separation_positions(&self, ub_rb_rows_numbers: &mut Vec<u32>, separation_positions: &mut Vec<Position>)
    {
        for bc in &self.boundary_conditions
        {
            if bc.is_type_same(BCType::Displacement)
            {
                let mut row = 0;
                for dof_parameter_data in self.state.nodes_dof_parameters_global.iter()
                {
                    if bc.is_dof_parameter_data_same(
                        dof_parameter_data.copy_dof_parameter(),         
                        dof_parameter_data.copy_node_number(),
                    )
                    {
                        separation_positions.push(Position(row, row));
                        ub_rb_rows_numbers.push(row);
                    }
                    row += 1;
                }
            }
        }
    }


    fn compose_ua_ra_rows_numbers(&self, ub_rb_rows_numbers: &[u32], ua_ra_rows_numbers: &mut Vec<u32>)
    {
        let mut i = 0;
        (0..self.state.nodes_dof_parameters_global.len()).for_each(|_|
            {
                if ub_rb_rows_numbers.iter().position(|n| *n == i).is_none()
                {
                    ua_ra_rows_numbers.push(i);
                }
                i += 1;
            });
    }


    fn compose_matrix_by_rows_numbers(&self, rows_numbers: &[u32], bc_type: BCType) -> Result<Vector<V>, String>
    {
        let mut all_elements = Vec::new();
        for row_number in rows_numbers
        {
            let converted_row_number = row_number as usize;

            let node_dof_parameter = self.state.nodes_dof_parameters_global[converted_row_number];
            if let Some(position) = self.boundary_conditions
                .iter()
                .position(|bc| bc
                    .is_dof_parameter_data_same(
                        node_dof_parameter.copy_dof_parameter(),
                        node_dof_parameter.copy_node_number(),
                    ) &&
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

        let mut converted_rows_numbers = rows_numbers.len();

        let matrix = Vector::create(&all_elements);
        Ok(matrix)
    }


    fn compose_displacements_matrix(
        &self, 
        ua_matrix: Matrix<V>,
        ub_matrix: Matrix<V>, 
        ua_ra_rows_numbers: &[u32],
        ub_rb_rows_numbers: &[u32]
    ) 
        -> Result<Matrix<V>, String>
    {
        let mut all_displacements_values = vec![V::from(0f32); self.state.nodes_dof_parameters_global.len()];

        let mut i = 0;
        let mut index = 0usize;
        while index < ua_ra_rows_numbers.len()
        {
            let displacement_value =  ua_matrix.get_element_value(&Position(i, 0))?;
            let converted_index = ua_ra_rows_numbers[index] as usize;
            all_displacements_values[converted_index] = displacement_value;
            i += 1;
            index += 1usize;
        }

        let mut j = 0;
        let mut index = 0usize;
        while index < ub_rb_rows_numbers.len()
        {
            let displacement_value = ub_matrix.get_element_value(&Position(j, 0))?;
            let converted_index = ub_rb_rows_numbers[index] as usize;
            all_displacements_values[converted_index] = displacement_value;
            j += 1;
            index += 1usize;
        }

        let mut rows_number = self.state.nodes_dof_parameters_global.len();

        let displacement_matrix = Matrix::create(
            rows_number, 1, &all_displacements_values,
        );
        Ok(displacement_matrix)
    }


    pub fn global_analysis(&mut self, stiffness_groups_update: bool) -> Result<GlobalAnalysisResult<V>, String>
    {
        if stiffness_groups_update
        {
            self.update_stiffness_groups()?;
        }
        self.update_nodes_dof_parameters_global()?;

        if self.boundary_conditions
            .iter()
            .position(|bc| bc.is_type_same(BCType::Displacement))
            .is_none()
        {
            return Err("FEModel: Model could not be analyzed because there are no restraints were \
                applied!".into())
        }

        let (mut global_stiffness_matrix, zero_rows_numbers, zero_columns_numbers) =
            self.compose_global_stiffness_matrix()?;

        let mut removed_zeros_rows_columns = Vec::new();
        for (zero_row_number, zero_column_number) in zero_rows_numbers
            .iter()
            .rev()
            .zip(zero_columns_numbers.iter().rev())
        {
            global_stiffness_matrix.remove_row(*zero_row_number);
            global_stiffness_matrix.remove_column(*zero_column_number);
            let matrix_element_position = Position(*zero_row_number, *zero_column_number);
            removed_zeros_rows_columns.push(matrix_element_position);
        }

        // let removed_zeros_rows_columns =
        //     global_stiffness_matrix.remove_zeros_rows_columns();

        self.shrink_of_nodes_dof_parameters(&removed_zeros_rows_columns)?;

        let mut ub_rb_rows_numbers = Vec::new();
        let mut separation_positions = Vec::new();
        self.compose_separation_positions(&mut ub_rb_rows_numbers, &mut separation_positions);

        let mut ua_ra_rows_numbers = Vec::new();
        self.compose_ua_ra_rows_numbers(&ub_rb_rows_numbers, &mut ua_ra_rows_numbers);

        let ra_matrix = self.compose_matrix_by_rows_numbers(&ua_ra_rows_numbers, BCType::Force)?;
        let ub_matrix = self.compose_matrix_by_rows_numbers(
            &ub_rb_rows_numbers, BCType::Displacement,
        )?;
        let rb_c_matrix = self.compose_matrix_by_rows_numbers(&ub_rb_rows_numbers, BCType::Force)?;

        let separated_matrix = separate(global_stiffness_matrix, separation_positions)?;

        let mut lhs_matrix = separated_matrix.ref_k_aa().clone();

        // lhs_matrix.try_to_symmetrize(self.state.tolerance);

        // let ua_matrix = lhs_matrix
        // .direct_solution(&ra_matrix.subtract_matrix(
        //     &separated_matrix.ref_k_ab().multiply_by_matrix(&ub_matrix)?)?, colsol_usage)?;

        let mut b = ra_matrix.subtract(&separated_matrix.ref_k_ab().multiply(&ub_matrix)?)?;
        let shape = b.get_shape();
        let mut v = Vec::new();
        for row in 0..shape.0
        {
            for column in 0..shape.1
            {
                let value = b.get_element_value(&Position(row, column))?;
                v.push(*value);
            }
        }
        let nn = v.len();
        let (a, maxa) = lhs_matrix.try_into_symmetric_compacted_matrix(self.state.rel_tol)?;
        factorization(&mut a, nn, &maxa)?;
        find_unknown(&a, &mut v, nn, &maxa);
        let ua_matrix = Vector::create(&v);

        let reactions_values_matrix = separated_matrix
            .ref_k_ba()
            .multiply(&ua_matrix)?
            .subtract(&separated_matrix.ref_k_bb().multiply(&ub_matrix)?)?
            .subtract(&rb_c_matrix)?;

        let reactions_values_matrix_shape = reactions_values_matrix.get_shape();
        let mut reactions_values = Vec::new();

        for row in 0..reactions_values_matrix_shape.0
        {
            for column in 0..reactions_values_matrix_shape.1
            {
                let reaction_value = reactions_values_matrix.get_element_value(&Position(row, column))?;
                reactions_values.push(*reaction_value);
            }
        }

        let mut reactions_dof_parameters_data = Vec::new();
        for row_number in ub_rb_rows_numbers.iter()
        {
            let converted_row_number = *row_number as usize;
            reactions_dof_parameters_data.push(self.state.nodes_dof_parameters_global[converted_row_number]);
        }
        let displacements_dof_parameters_data = self.state.nodes_dof_parameters_global.clone();
        let displacements_values_matrix = self.compose_displacements_matrix(
            ua_matrix, ub_matrix, &ua_ra_rows_numbers, &ub_rb_rows_numbers,
        )?;

        let displacements_values_matrix_shape = displacements_values_matrix.get_shape();
        let mut displacements_values = Vec::new();

        for row in 0..displacements_values_matrix_shape.0
        {
            for column in 0..displacements_values_matrix_shape.1
            {
                let displacement_value = displacements_values_matrix
                    .get_element_value(&Position(row, column))?;
                displacements_values.push(*displacement_value);
            }
        }

        let global_analysis_result = GlobalAnalysisResult::create(
            reactions_values, reactions_dof_parameters_data,  displacements_values, displacements_dof_parameters_data,
        );

        Ok(global_analysis_result)
    }


    pub fn prepare_to_global_analysis(&mut self, stiffness_groups_update: bool) -> Result<(), String>
    {
        if stiffness_groups_update
        {
            self.update_stiffness_groups()?;
        }
        self.update_nodes_dof_parameters_global()?;

        if self.boundary_conditions
            .iter()
            .position(|bc| bc.is_type_same(BCType::Displacement))
            .is_none()
        {
            return Err("FEModel: Model could not be analyzed because there are no restraints were applied!".into())
        }

        let (mut global_stiffness_matrix, zero_rows_numbers, zero_columns_numbers) =
            self.compose_global_stiffness_matrix()?;

        let mut removed_zeros_rows_columns = Vec::new();
        for (zero_row_number, zero_column_number) in zero_rows_numbers
            .iter()
            .rev()
            .zip(zero_columns_numbers.iter().rev())
        {
            global_stiffness_matrix.remove_selected_row(*zero_row_number);
            global_stiffness_matrix.remove_selected_column(*zero_column_number);
            let matrix_element_position = Position(*zero_row_number, *zero_column_number);
            removed_zeros_rows_columns.push(matrix_element_position);
        }

        // let removed_zeros_rows_columns =
        //     global_stiffness_matrix.remove_zeros_rows_columns();

        self.shrink_of_nodes_dof_parameters(&removed_zeros_rows_columns)?;

        let mut ub_rb_rows_numbers = Vec::new();
        let mut separation_positions = Vec::new();
        self.compose_separation_positions(&mut ub_rb_rows_numbers, &mut separation_positions);

        let mut ua_ra_rows_numbers = Vec::new();
        self.compose_ua_ra_rows_numbers(&ub_rb_rows_numbers, &mut ua_ra_rows_numbers);

        let ra_matrix = self.compose_matrix_by_rows_numbers(&ua_ra_rows_numbers, BCType::Force)?;
        let ub_matrix = self.compose_matrix_by_rows_numbers(
            &ub_rb_rows_numbers, BCType::Displacement,
        )?;
        let rb_c_matrix = self.compose_matrix_by_rows_numbers(&ub_rb_rows_numbers, BCType::Force)?;

        let separated_matrix = separate(global_stiffness_matrix, separation_positions)?;

        self.state.optional_ua_ra_rows_numbers = Some(ua_ra_rows_numbers);
        self.state.optional_ub_rb_rows_numbers = Some(ub_rb_rows_numbers);
        self.state.optional_ub_matrix = Some(ub_matrix);
        self.state.optional_ra_matrix = Some(ra_matrix);
        self.state.optional_rb_c_matrix = Some(rb_c_matrix);
        self.state.optional_separated_matrix = Some(separated_matrix);

        Ok(())
    }


    pub fn calculate_ua_matrix(&mut self, colsol_usage: bool, f: fn(data: &str)) -> Result<(), String>
    {
        if self.state.optional_ra_matrix.is_some() && 
            self.state.optional_ub_matrix.is_some() &&
            self.state.optional_separated_matrix.is_some()
        {
            let mut lhs_matrix = self.state.optional_separated_matrix
                .as_ref()
                .unwrap()
                .ref_k_aa()
                .clone();
            let rhs_matrix = self.state.optional_ra_matrix
                .as_ref()
                .unwrap()
                .subtract(&self.state.optional_separated_matrix
                    .as_ref()
                    .unwrap()
                    .ref_k_ab()
                    .multiply(&self.state.optional_ub_matrix.as_ref().unwrap())?)?;

            // lhs_matrix.try_to_symmetrize(self.state.tolerance);
            // f(&format!("{:?}, colsol usage: {colsol_usage}", lhs_matrix.ref_matrix_type()));
            // let ua_matrix = lhs_matrix.direct_solution(&rhs_matrix, colsol_usage)?;

            let shape = rhs_matrix.get_shape();
            let mut v = Vec::new();
            for row in 0..shape.0
            {
                for column in 0..shape.1
                {
                    let value = rhs_matrix.get_element_value(&Position(row, column))?;
                    v.push(*value);
                }
            }
            let nn = v.len();
            let (a, maxa) = lhs_matrix.try_into_symmetric_compacted_matrix(self.state.rel_tol)?;
            factorization(&mut a, nn, &maxa)?;
            find_unknown(&a, &mut v, nn, &maxa);
            let ua_matrix = Vector::create(&v);

            self.state.optional_ua_matrix = Some(ua_matrix);
            Ok(())
        }
        else
        {
            let error_message = "FEModel: Calculation of Ua matrix: \
                Insufficient data to calculate Ua matrix!";
            Err(error_message.to_string())
        }
    }


    pub fn extract_data_for_ua_matrix_calculation(&mut self) -> Result<(usize, usize, Vec<V>, Vec<V>), String>
    {
        if self.state.optional_ra_matrix.is_some() && 
            self.state.optional_ub_matrix.is_some() &&
            self.state.optional_separated_matrix.is_some()
        {
            let lhs_matrix = self.state.optional_separated_matrix.as_ref().unwrap().ref_k_aa();
            let rhs_matrix = self.state.optional_ra_matrix
                .as_ref()
                .unwrap()
                .subtract(&self.state.optional_separated_matrix
                    .as_ref()
                    .unwrap()
                    .ref_k_ab()
                    .multiply(&self.state.optional_ub_matrix.as_ref().unwrap())?)?;
            
            let lhs_matrix_shape = lhs_matrix.copy_shape();
            let a_rows_number = lhs_matrix_shape.0;
            let a_columns_number = lhs_matrix_shape.1;

            let mut a_elements_values = vec![V::from(0f32); a_rows_number * a_columns_number];

            for (element_position, element_value) in lhs_matrix.get_elements()
            {
                let row_position = element_position.0;
                let column_position = element_position.1;
                let index = row_position * a_columns_number + column_position;
                a_elements_values[index] = *element_value;
                if row_position != column_position
                {
                    let symmetric_index = column_position * a_columns_number + row_position;
                    a_elements_values[symmetric_index] = *element_value;
                }
            }


            let mut b_elements_values = vec![V::from(0f32); a_rows_number];

            for (element_position, element_value) in rhs_matrix.get_elements()
            {
                let row_position = element_position.0;
                let column_position = element_position.1;
                let index = row_position + column_position;
                b_elements_values[index] = *element_value;
            }

            Ok((a_rows_number, a_columns_number, a_elements_values, b_elements_values))
        }
        else
        {
            Err("FEModel: Calculation of Ua matrix: Insufficient data to calculate Ua matrix!".to_string())
        }
    }


    pub fn receive_ua_matrix_data(
        &mut self, 
        a_rows_number: usize, 
        a_columns_number: usize, 
        a_elements_values: Vec<V>,
    ) 
        -> Result<(), String>
    {
        if a_elements_values.len() != a_rows_number * a_columns_number
        {
            Err("FEModel: Receiving of Ua matrix: Inapropriate number of elements!".to_string()) 
        }
        let ua_matrix = Vector::create(&a_elements_values);
        self.state.optional_ua_matrix = Some(ua_matrix);
        Ok(())
    }


    pub fn extract_global_analysis_result(&mut self) -> Result<GlobalAnalysisResult<V>, String>
    {
        if self.state.optional_ua_matrix.is_none() || 
            self.state.optional_ub_matrix.is_none() ||
            self.state.optional_separated_matrix.is_none() ||
            self.state.optional_rb_c_matrix.is_none() ||
            self.state.optional_ua_ra_rows_numbers.is_none() ||
            self.state.optional_ub_rb_rows_numbers.is_none()
        {
            let error_message = "FEModel: Extraction of global analysis data: \
                Insufficient data to perform global analysis!";
            return Err(error_message.to_string());
        }

        let separated_matrix = self.state.optional_separated_matrix.as_ref().unwrap();
        let ua_matrix = self.state.optional_ua_matrix.as_ref().unwrap();
        let ub_matrix = self.state.optional_ub_matrix.as_ref().unwrap();
        let rb_c_matrix = self.state.optional_rb_c_matrix.as_ref().unwrap();
        let ua_ra_rows_numbers = self.state.optional_ua_ra_rows_numbers.as_ref().unwrap();
        let ub_rb_rows_numbers = self.state.optional_ub_rb_rows_numbers.as_ref().unwrap();

        let reactions_values_matrix = separated_matrix
            .ref_k_ba()
            .multiply(&ua_matrix)?
            .add(&separated_matrix.ref_k_bb().multiply(&ub_matrix)?)?
            .subtract(&rb_c_matrix)?;

        let reactions_values_matrix_shape = reactions_values_matrix.get_shape();
        let mut reactions_values = Vec::new();

        for row in 0..reactions_values_matrix_shape.0
        {
            for column in 0..reactions_values_matrix_shape.1
            {
                let reaction_value = reactions_values_matrix.get_element_value(&Position(row, column))?;
                reactions_values.push(*reaction_value);
            }
        }

        let mut reactions_dof_parameters_data = Vec::new();
        for row_number in ub_rb_rows_numbers
        {
            reactions_dof_parameters_data.push(self.state.nodes_dof_parameters_global[row_number as usize]);
        }
        let displacements_dof_parameters_data = self.state.nodes_dof_parameters_global.clone();
        let displacements_values_matrix = self.compose_displacements_matrix(
            ua_matrix.clone(), ub_matrix.clone(), &ua_ra_rows_numbers, &ub_rb_rows_numbers,
        )?;

        let displacements_values_matrix_shape = displacements_values_matrix.get_shape();
        let mut displacements_values = Vec::new();

        for row in 0..displacements_values_matrix_shape.0
        {
            for column in 0..displacements_values_matrix_shape.1
            {
                let displacement_value = displacements_values_matrix
                    .get_element_value(&Position(row, column))?;
                displacements_values.push(*displacement_value);
            }
        }

        let global_analysis_result = GlobalAnalysisResult::create(
            reactions_values, reactions_dof_parameters_data, displacements_values, displacements_dof_parameters_data,
        );

        self.reset_optional_state_values();

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


    pub fn extract_unique_elements_of_rotation_matrix(&self, element_number: &u32)
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
        if let Some(position) = self.boundary_conditions
            .iter()
            .position(|bc| bc.is_type_same(bc_type) && bc.is_number_same(number))
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


    pub fn convert_uniformly_distributed_surface_force_to_nodal_forces(&self, element_number: T,
        uniformly_distributed_surface_force_value: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        if let Some(finite_element) = self.elements.get(&element_number) 
        {
            match finite_element.copy_fe_type()
            {
                FEType::Plate4n4ip => 
                {
                    if let Some(element) = finite_element.element.as_any().downcast_ref::<Plate4n4ip<T, V>>()
                    {
                        element.convert_uniformly_distributed_surface_force_to_nodal_forces(uniformly_distributed_surface_force_value, 
                            &self.nodes, self.state.tolerance)
                    }
                    else 
                    {
                        return Err(format!("FEModel: Incorrect type of element {:?}!", element_number));
                    }
                },
                _ => 
                {
                    return Err(format!("FEModel: Incorrect type of element {:?}!", element_number));
                },
            }
        }
        else
        {
            Err(format!("FEModel: Finite element with number {:?} does not exist!", element_number))
        }
    }


    pub fn convert_uniformly_distributed_line_force_to_nodal_forces(&self, element_number: T,
        uniformly_distributed_line_force_value: V) -> Result<ExtendedMatrix<T, V>, String>
    {
        if let Some(finite_element) = self.elements.get(&element_number) 
        {
            match finite_element.copy_fe_type()
            {
                FEType::Beam2n1ipT => 
                {
                    if let Some(element) = finite_element.element.as_any().downcast_ref::<Beam2n1ipT<T, V>>()
                    {
                        element.convert_uniformly_distributed_line_force_to_nodal_forces(uniformly_distributed_line_force_value, 
                            &self.nodes, self.state.tolerance)
                    }
                    else 
                    {
                        return Err(format!("FEModel: Incorrect type of element {:?}!", element_number));
                    }
                },
                _ => 
                {
                    return Err(format!("FEModel: Incorrect type of element {:?}!", element_number));
                },
            }
        }
        else
        {
            Err(format!("FEModel: Finite element with number {:?} does not exist!", element_number))
        }
    }
}
