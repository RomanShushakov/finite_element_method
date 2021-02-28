use crate::fem::
    {
        FeNode, FEData, FiniteElement, StiffnessGroup, DOFParameterData, BoundaryCondition,
        GlobalAnalysisResult
    };
use crate::fem::{FEType, GlobalDOFParameter, BCType};
use crate::fem::compose_stiffness_sub_groups;
use crate::fem::GLOBAL_DOF;

use crate::extended_matrix::{ExtendedMatrix, MatrixElementPosition, ZerosRowColumn};
use crate::extended_matrix::Operation;
use crate::extended_matrix::extract_element_value;

use crate::{ElementsNumbers, ElementsValues};

use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::rc::Rc;
use std::cell::RefCell;
use std::collections::HashSet;
use std::iter::FromIterator;


pub struct SeparatedMatrix<T, V>
{
    pub k_aa: ExtendedMatrix<T, V>,
    pub k_ab: ExtendedMatrix<T, V>,
    pub k_ba: ExtendedMatrix<T, V>,
    pub k_bb: ExtendedMatrix<T, V>,
}


pub struct State<T>
{
    pub stiffness_groups: Vec<StiffnessGroup<T>>,
    pub nodes_dof_parameters_global: Vec<DOFParameterData<T>>,
}


pub struct FEModel<T, V>
{
    pub nodes: Vec<Rc<RefCell<FeNode<T, V>>>>,
    pub elements: Vec<FiniteElement<T, V>>,
    pub boundary_conditions: Vec<BoundaryCondition<T, V>>,
    pub state: State<T>,
}


impl<T, V> FEModel<T, V>
    where T: Copy + PartialEq + Into<ElementsNumbers> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + From<ElementsNumbers> + Eq + Hash + SubAssign + Debug +
             Mul<Output = T> + PartialOrd + Default + Add<Output = T> + AddAssign + 'static,
          V: Copy + From<ElementsValues> + Sub<Output = V> + Default + Mul<Output = V> +
             Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign +
             SubAssign + Into<ElementsValues> + 'static
{
    pub fn create() -> Self
    {
        let state = State { stiffness_groups: Vec::new(), nodes_dof_parameters_global: Vec::new() };
        FEModel { nodes: Vec::new(), elements: Vec::new(), boundary_conditions: Vec::new(), state }
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
            for node in self.nodes.iter()
            {
                nodes_numbers.push(node.borrow().number);
            }
            let mut position = T::default();
            let columns_number = T::from(nodes_numbers.len() as ElementsNumbers);
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
                        position += T::from(1);
                    }
                    let stiffness_sub_groups =
                         compose_stiffness_sub_groups(position,
                        columns_number, excluded,
                        v_lhs[j])?;
                    stiffness_groups.extend(stiffness_sub_groups);
                    position += T::from(1);
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
                position += T::from(1);
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
        for node in &self.nodes
        {
            for dof in 0..GLOBAL_DOF
            {
                let dof_parameter =
                    GlobalDOFParameter::iterator().nth(dof as usize)
                        .ok_or("FEModel: Could not find dof parameter!")?;
                let dof_parameter_data = DOFParameterData {
                    node_number: node.as_ref().borrow().number,
                    dof_parameter: *dof_parameter
                };
                nodes_dof_parameters.push(dof_parameter_data);
            }
        }
        self.state.nodes_dof_parameters_global = nodes_dof_parameters;
        Ok(())
    }


    pub fn add_node(&mut self, number: T, x: V, y: V, z: V) -> Result<(), String>
    {
        if self.nodes.iter().position(|node|
            node.as_ref().borrow().number_same(number) ||
            node.as_ref().borrow().coordinates_same(x, y, z)).is_none()
        {
            let node = FeNode::create(number, x, y, z);
            self.nodes.push(Rc::new(RefCell::new(node)));
            self.update_stiffness_groups()?;
            Ok(())
        }
        else
        {
            Err(format!("FEModel: Node {} could not be added because node with same number or \
                with the same coordinates does already exist!", number.into()))
        }
    }


    pub fn update_node(&mut self, number: T, x: V, y: V, z: V) -> Result<(), String>
    {
        if self.nodes.iter().position(|node|
            !node.as_ref().borrow().number_same(number) &&
            node.as_ref().borrow().coordinates_same(x, y, z)).is_some()
        {
            return Err(format!("FEModel: Node {} could not be updated because the node with the \
                same coordinates does already exist!", number.into()))
        }
        if let Some(position) = self.nodes.iter().position(|node|
            node.as_ref().borrow().number_same(number))
        {
            self.nodes[position].borrow_mut().update(x, y, z);
            for element in self.elements
                .iter_mut()
                .filter(|element|
                    element.node_belong_element(number))
            {
                element.refresh()?;
            }
            return Ok(());
        }
        Err(format!("FEModel: Node {} could not be updated because it does not \
            exist!", number.into()))
    }


    pub fn delete_node(&mut self, number: T) -> Result<(), String>
    {
        if let Some(position) = self.nodes.iter().position(|node|
            node.as_ref().borrow().number_same(number))
        {
            while let Some(position) = self.elements
                .iter()
                .position(|element|
                    element.node_belong_element(number))
            {
                self.elements.remove(position);
            }
            while let Some(position) = self.boundary_conditions
                .iter()
                .position(|bc|
                    bc.node_number_same(number))
            {
                self.boundary_conditions.remove(position);
            }
            self.nodes.remove(position);
            self.update_stiffness_groups()?;
            return Ok(());
        }
        Err(format!("FEModel: Node {} could not be deleted because it does not exist!",
                    number.into()))
    }


    pub fn add_element(&mut self, element_type: FEType, nodes_numbers: Vec<T>,
        mut data: FEData<T, V>) -> Result<(), String>
    {
        if self.elements.iter().position(|element|
            element.number_same(data.number)).is_some()
        {
            return Err(format!("FEModel: Element {} could not be added! The element with the same \
             number does already exist!", data.number.into()));
        }
        let nodes_numbers_set = HashSet::<T>::from_iter(
            nodes_numbers.iter().cloned());
        if nodes_numbers.len() != nodes_numbers_set.len()
        {
            return Err(format!("FEModel: Element {} could not be added! All nodes numbers \
                should be unique!", data.number.into()));
        }
        if self.elements.iter().position(|element|
            element.element_type == element_type &&
            element.nodes_numbers_same(nodes_numbers.clone())).is_some()
        {
            return Err(format!("FEModel: Element {} could not be added! The element with the same \
                type and with same nodes numbers does already exist!", data.number.into()));
        }
        for node_number in nodes_numbers.iter()
        {
            if let Some(position) = self.nodes.iter().position(|node|
                node.as_ref().borrow().number_same(*node_number))
            {
                data.nodes.push(Rc::clone(&self.nodes[position]));
            }
        }
        if nodes_numbers.len() == data.nodes.len()
        {
            let element = FiniteElement::create(element_type, data)?;
            self.elements.push(element);
        }
        else
        {
            return Err(format!("FEModel: Element {} could not be added! Some node does not exist!",
                               data.number.into()));
        }
        Ok(())
    }


    pub fn update_element(&mut self, nodes_numbers: Vec<T>,
        mut data: FEData<T, V>) -> Result<(), String>
    {
        let nodes_numbers_set = HashSet::<T>::from_iter(
            nodes_numbers.iter().cloned());
        if nodes_numbers.len() != nodes_numbers_set.len()
        {
            return Err(format!("FEModel: Element {} could not be updated! All nodes numbers \
                should be unique!", data.number.into()));
        }
        for node_number in nodes_numbers.iter()
        {
            if let Some(position) = self.nodes.iter().position(|node|
                node.as_ref().borrow().number_same(*node_number))
            {
                data.nodes.push(Rc::clone(&self.nodes[position]));
            }
        }
        if self.elements.iter().position(|element|
            !element.number_same(data.number) &&
            element.nodes_numbers_same(nodes_numbers.clone())).is_some()
        {
            return Err(format!("FEModel: Element {} could not be added! The element with the same \
                nodes numbers does already exist!", data.number.into()));
        }
        if nodes_numbers.len() == data.nodes.len()
        {
            if let Some(position) = self.elements.iter().position(|element|
                element.number_same(data.number))
            {
                self.elements[position].update(data)?;
            }
            else
            {
               return Err(format!("FEModel: Element {} could not be updated because it does not \
                exist!", data.number.into()));
            }

        }
        else
        {
            return Err(format!("FEModel: Element {} could not be updated! Some node does not exist!",
                               data.number.into()));
        }
        Ok(())
    }


    pub fn delete_element(&mut self, number: T) -> Result<(), String>
    {
        if let Some(position) = self.elements.iter().position(|element|
            element.number_same(number))
        {
            self.elements.remove(position);
            return Ok(());
        }
        Err(format!("FEModel: Element {} could not be deleted because it does not exist!",
                    number.into()))
    }


    fn compose_global_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
    {
        if self.elements.is_empty()
        {
            return Err("FEModel: Global stiffness matrix could not be composed because there are \
                no elements in the model!");
        }
        if self.nodes.iter().any(|node|
            self.elements.iter().position(|element|
                element.node_belong_element(node.as_ref().borrow().number)).is_none())
        {
            return Err("FEModel: Global stiffness matrix could not be composed because there are \
                free nodes does exist!");
        }
        let mut global_stiffness_matrix = ExtendedMatrix::create(
            T::from(self.nodes.len() as ElementsNumbers * GLOBAL_DOF),
            T::from(self.nodes.len() as ElementsNumbers * GLOBAL_DOF),
            vec![V::default(); (self.nodes.len() as ElementsNumbers * GLOBAL_DOF *
                self.nodes.len() as ElementsNumbers * GLOBAL_DOF) as usize]);
        for element in &self.elements
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
            bc.number_same(number) && bc.type_same(bc_type)).is_some()
        {
            return Err(format!("FEModel: {} could not be added because the same {} number does \
                already exist!", bc_type.as_str(), bc_type.as_str().to_lowercase()));
        }
        if self.boundary_conditions.iter().position(|bc|
            bc.dof_parameter_data_same(dof_parameter, node_number)).is_some()
        {
            return Err(format!("FEModel: {} could not be added because the the force or \
                displacement with the same dof parameter data does already exist!",
                               bc_type.as_str()));
        }
        if self.nodes.iter().position(|node|
            node.as_ref().borrow().number_same(node_number)).is_none()
        {
            return Err(format!("FEModel: {} could not be added because the current node number \
                does not exist!", bc_type.as_str()));
        }
        let bc = BoundaryCondition::create(
            bc_type, number, node_number, dof_parameter, value);
        self.boundary_conditions.push(bc);
        Ok(())
    }


    pub fn update_bc(&mut self, bc_type: BCType, number: T, node_number: T,
        dof_parameter: GlobalDOFParameter, value: V) -> Result<(), String>
    {
        if self.nodes.iter().position(|node|
            node.as_ref().borrow().number_same(node_number)).is_none()
        {
            return Err(format!("FEModel: {} could not be updated because the current node number \
                does not exist!", bc_type.as_str()));
        }
        if self.boundary_conditions.iter().position(|bc|
                bc.dof_parameter_data_same(
                    dof_parameter, node_number) && !bc.number_same(number)).is_some()
        {
            return Err(format!("FEModel: {} could not be updated because the the force or \
                displacement with the same dof parameter data does already exist!",
                               bc_type.as_str()));
        }
        if let Some(position) =  self.boundary_conditions.iter().position(|bc|
            bc.number_same(number) && bc.type_same(bc_type))
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
            bc.number_same(number) && bc.type_same(bc_type))
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
            let dof_parameter_data =
                self.state.nodes_dof_parameters_global.remove(
                    row_column.column.into() as usize);
            if let Some(position) = self.boundary_conditions
                .iter()
                .position(|bc|
                    bc.dof_parameter_data_same(
                        dof_parameter_data.dof_parameter, dof_parameter_data.node_number))
            {
                let bc_type = self.boundary_conditions[position].extract_bc_type();
                let dof_parameter = dof_parameter_data.dof_parameter;
                let node_number = dof_parameter_data.node_number;
                return Err(format!("FEModel: Model could not be analyzed because where are \
                    no stiffness to withstand {}::{:?} applied at node {:?}!",
                                   bc_type.as_str(), dof_parameter, node_number))
            }
        }
        Ok(())
    }


    fn compose_separation_positions(&self, ub_rb_rows_numbers: &mut Vec<T>,
        separation_positions: &mut Vec<MatrixElementPosition<T>>)
    {
        for bc in &self.boundary_conditions
        {
            if bc.type_same(BCType::Displacement)
            {
                let mut row = T::default();
                for dof_parameter_data in
                    &self.state.nodes_dof_parameters_global
                {
                    if bc.dof_parameter_data_same(
                        dof_parameter_data.dof_parameter, dof_parameter_data.node_number)
                    {
                        separation_positions.push(MatrixElementPosition { row, column: row });
                        ub_rb_rows_numbers.push(row);
                    }
                    row += T::from(1);
                }
            }
        }
    }


    fn compose_ua_ra_rows_numbers(&self, ub_rb_rows_numbers: &Vec<T>,
        ua_ra_rows_numbers: &mut Vec<T>)
    {
        for i in 0..self.state.nodes_dof_parameters_global.len()
        {
            if ub_rb_rows_numbers.iter().position(|n|
                *n == T::from(i as ElementsNumbers)).is_none()
            {
                ua_ra_rows_numbers.push(T::from(i as ElementsNumbers));
            }
        }
    }


    fn compose_matrix_by_rows_numbers(&self, rows_numbers: &Vec<T>) -> ExtendedMatrix<T, V>
    {
        let mut all_elements = Vec::new();
        for row_number in rows_numbers
        {
            let node_dof_parameter =
                self.state.nodes_dof_parameters_global[(*row_number).into() as usize];
            if let Some(position) = self.boundary_conditions
                .iter()
                .position(|bc|
                    bc.dof_parameter_data_same(
                        node_dof_parameter.dof_parameter, node_dof_parameter.node_number))
            {
                let value = self.boundary_conditions[position].extract_value();
                all_elements.push(value);
            }
            else
            {
                all_elements.push(V::default());
            }
        }
        let matrix = ExtendedMatrix::create(
            T::from(rows_numbers.len() as ElementsNumbers),
            T::from(1),
            all_elements);
        matrix
    }


    fn compose_displacements_matrix(&self, ua_matrix: ExtendedMatrix<T, V>,
        ub_matrix: ExtendedMatrix<T, V>, ua_ra_rows_numbers: &Vec<T>,
        ub_rb_rows_numbers: &Vec<T>) -> ExtendedMatrix<T, V>
    {
        let ua_values = ua_matrix.extract_all_elements_values();
        let ub_values = ub_matrix.extract_all_elements_values();
        let mut all_displacements_values =
            vec![V::default(); self.state.nodes_dof_parameters_global.len()];
        for i in  0..ua_ra_rows_numbers.len()
        {
            let displacement_value = extract_element_value(
                T::from(i as ElementsNumbers), T::from(0), &ua_values);
            all_displacements_values[ua_ra_rows_numbers[i].into() as usize] = displacement_value;
        }
        for j in  0..ub_values.len()
        {
            let displacement_value = extract_element_value(
                T::from(j as ElementsNumbers), T::from(0), &ub_values);
            all_displacements_values[ub_rb_rows_numbers[j].into() as usize] = displacement_value;
        }
        let rows_number =
            T::from(self.state.nodes_dof_parameters_global.len() as ElementsNumbers);
        let displacement_matrix =
            ExtendedMatrix::create(
                rows_number, T::from(1), all_displacements_values);
        displacement_matrix
    }


    pub fn global_analysis(&mut self) -> Result<GlobalAnalysisResult<T, V>, String>
    {
        self.update_nodes_dof_parameters_global()?;
        if self.boundary_conditions.iter().position(|bc|
            bc.type_same(BCType::Displacement)).is_none()
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
            global_stiffness_matrix.separate(separation_positions)?;
        let ua_matrix = separated_matrix.k_aa
            .naive_gauss_elimination(&ra_matrix.add_subtract_matrix(
            &separated_matrix.k_ab.multiply_by_matrix(&ub_matrix)?,
            Operation::Subtraction)?)?;
        let reactions_values_matrix = separated_matrix.k_ba
            .multiply_by_matrix(&ua_matrix)?
            .add_subtract_matrix(
                &separated_matrix.k_bb
                    .multiply_by_matrix(&ub_matrix)?, Operation::Addition)?;
        let all_reactions =
            reactions_values_matrix.extract_all_elements_values();
        let reactions_values_matrix_shape = reactions_values_matrix.get_shape();
        let mut reactions_values = Vec::new();
        for row in 0..reactions_values_matrix_shape.0.into()
        {
            for column in 0..reactions_values_matrix_shape.1.into()
            {
                let reaction_value =
                    extract_element_value(T::from(row),
                        T::from(column), &all_reactions);
                reactions_values.push(reaction_value);
            }
        }
        let mut reactions_dof_parameters_data = Vec::new();
        for row_number in &ub_rb_rows_numbers
        {
            reactions_dof_parameters_data.push(
                self.state.nodes_dof_parameters_global[(*row_number).into() as usize]);
        }
        let displacements_dof_parameters_data =
            self.state.nodes_dof_parameters_global.clone();
        let displacements_values_matrix = self.compose_displacements_matrix(
            ua_matrix, ub_matrix, &ua_ra_rows_numbers, &ub_rb_rows_numbers);
        let all_displacements =
            displacements_values_matrix.extract_all_elements_values();
        let displacements_values_matrix_shape = displacements_values_matrix.get_shape();
        let mut displacements_values = Vec::new();
        for row in 0..displacements_values_matrix_shape.0.into()
        {
            for column in 0..displacements_values_matrix_shape.1.into()
            {
                let displacement_value =
                    extract_element_value(T::from(row),
                        T::from(column), &all_displacements);
                displacements_values.push(displacement_value);
            }
        }
        let global_analysis_result =
            GlobalAnalysisResult::create(
                reactions_values, reactions_dof_parameters_data,
                displacements_values, displacements_dof_parameters_data);
        Ok(global_analysis_result)
    }
}
