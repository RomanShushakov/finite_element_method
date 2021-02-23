use crate::fem::{FeNode, FEData, FiniteElement, StiffnessGroup, ForceBC, DisplacementBC};
use crate::fem::{FEType, GlobalForceDisplacementComponent};
use crate::fem::compose_stiffness_sub_groups;
use crate::{ElementsNumbers, ElementsValues};
use crate::extended_matrix::ExtendedMatrix;

use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::rc::Rc;
use std::cell::RefCell;
use std::collections::HashSet;
use std::iter::FromIterator;


pub const GLOBAL_DOF: ElementsNumbers = 6;


pub struct FEModel<T, V>
{
    pub nodes: Vec<Rc<RefCell<FeNode<T, V>>>>,
    pub elements: Vec<FiniteElement<T, V>>,
    pub stiffness_groups: Vec<StiffnessGroup<T>>,
    pub applied_loads: Vec<ForceBC<T, V>>,
    pub applied_displacements: Vec<DisplacementBC<T, V>>,
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
        FEModel { nodes: Vec::new(), elements: Vec::new(), stiffness_groups: Vec::new(),
            applied_loads: Vec::new(), applied_displacements: Vec::new() }
    }


    fn update_stiffness_groups(&mut self) -> Result<(), &str>
    {
        let mut stiffness_groups = Vec::new();
        if self.nodes.len() < 2
        {
            self.stiffness_groups = Vec::new();
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
        self.stiffness_groups = stiffness_groups;
        Ok(())
    }


    pub fn add_node(&mut self, number: T, x: V, y: V, z: V) -> Result<(), String>
    {
        if self.nodes.iter().position(|node|
            node.as_ref().borrow().number == number).is_none()
        {
            let node = FeNode::create(number, x, y, z);
            self.nodes.push(Rc::new(RefCell::new(node)));
            self.update_stiffness_groups()?;
            Ok(())
        }
        else
        {
            Err(format!("FEModel: Node {} could not be added because it does already exist!",
                        number.into()))
        }
    }


    pub fn update_node(&mut self, number: T, x: V, y: V, z: V) -> Result<(), String>
    {
        if let Some(position) = self.nodes.iter().position(|node|
            node.as_ref().borrow().number == number)
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
        Err(format!("FEModel: Node {} could not be updated!", number.into()))
    }


    pub fn delete_node(&mut self, number: T) -> Result<(), String>
    {
        if let Some(position) = self.nodes.iter().position(|node|
            node.as_ref().borrow().number == number)
        {
            while let Some(position) = self.elements
                .iter()
                .position(|element|
                    element.node_belong_element(number))
            {
                self.elements.remove(position);
            }
            while let Some(position) = self.applied_loads
                .iter()
                .position(|force_bc|
                    force_bc.force.node_number == number)
            {
                self.applied_loads.remove(position);
            }
            while let Some(position) = self.applied_displacements
                .iter()
                .position(|displacement_bc|
                    displacement_bc.displacement.node_number == number)
            {
                self.applied_displacements.remove(position);
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
                node.as_ref().borrow().number == *node_number)
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
                node.as_ref().borrow().number == *node_number)
            {
                data.nodes.push(Rc::clone(&self.nodes[position]));
            }
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


    pub fn compose_global_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
    {
        if self.elements.is_empty()
        {
            return Err("FEModel: Global stiffness matrix could not be composed because there are \
                no elements in the model!");
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
                if let Some(position) = self.stiffness_groups
                    .iter()
                    .position(|group|
                        { group.stiffness_type == element_stiffness_group.stiffness_type &&
                        group.number_1 == element_stiffness_group.number_1 &&
                        group.number_2 == element_stiffness_group.number_2 })
                {
                    global_stiffness_matrix.add_sub_matrix(
                        &element_stiffness_matrix,
                        &self.stiffness_groups[position].positions,
                        &element_stiffness_group.positions);
                }
            }
        }
        Ok(global_stiffness_matrix)
    }


    pub fn add_load(&mut self, number: T, node_number: T,
        component: GlobalForceDisplacementComponent, value: V) -> Result<(), &str>
    {
        if self.applied_loads.iter().position(|f|
            f.number == number).is_some()
        {
            return Err("FEModel: Force could not be added because the same force number does \
                already exist!");
        }
        if self.nodes.iter().position(|node|
            node.as_ref().borrow().number == node_number).is_none()
        {
            return Err("FEModel: Force could not be added because the current node number does \
                not exist!");
        }
        let force_bc = ForceBC::create(number, node_number, component, value);
        self.applied_loads.push(force_bc);
        Ok(())
    }


    pub fn update_load(&mut self, number: T, node_number: T,
        component: GlobalForceDisplacementComponent, value: V) -> Result<(), &str>
    {
        if self.nodes.iter().position(|node|
            node.as_ref().borrow().number == node_number).is_none()
        {
            return Err("FEModel: Force could not be updated because the current node number does \
                not exist!");
        }
        if let Some(position) =  self.applied_loads.iter().position(|f|
            f.number == number)
        {

            self.applied_loads[position].update(node_number, component, value);
            Ok(())
        }
        else
        {
            Err("FEModel: Force could not be updated because current force number does not exist!")
        }
    }


    pub fn delete_load(&mut self, number: T) -> Result<(), &str>
    {
        if let Some(position) =  self.applied_loads.iter().position(|f|
            f.number == number)
        {
            self.applied_loads.remove(position);
            Ok(())
        }
        else
        {
            Err("FEModel: Force could not be deleted because current force number does not exist!")
        }
    }


    pub fn add_displacement(&mut self, number: T, node_number: T,
        component: GlobalForceDisplacementComponent, value: V) -> Result<(), &str>
    {
        if self.applied_displacements.iter().position(|d|
            d.number == number).is_some()
        {
            return Err("FEModel: Displacement could not be added because the same displacement \
                number does already exist!");
        }
        if self.nodes.iter().position(|node|
            node.as_ref().borrow().number == node_number).is_none()
        {
            return Err("FEModel: Displacement could not be added because the current node number \
                does not exist!");
        }
        let displacement_bc = DisplacementBC::create(
            number, node_number, component, value);
        self.applied_displacements.push(displacement_bc);
        Ok(())
    }


    pub fn update_displacement(&mut self, number: T, node_number: T,
        component: GlobalForceDisplacementComponent, value: V) -> Result<(), &str>
    {
        if self.nodes.iter().position(|node|
            node.as_ref().borrow().number == node_number).is_none()
        {
            return Err("FEModel: Displacement could not be updated because the current node number \
                does not exist!");
        }
        if let Some(position) =  self.applied_displacements.iter().position(|d|
            d.number == number)
        {
            self.applied_displacements[position].update(node_number, component, value);
            Ok(())
        }
        else
        {
            Err("FEModel: Displacement could not be updated because current displacement number \
                does not exist!")
        }
    }


    pub fn delete_displacement(&mut self, number: T) -> Result<(), &str>
    {
        if let Some(position) =  self.applied_displacements.iter().position(|d|
            d.number == number)
        {
            self.applied_displacements.remove(position);
            Ok(())
        }
        else
        {
            Err("FEModel: Displacement could not be deleted because current displacement number \
                does not exist!")
        }
    }
}
