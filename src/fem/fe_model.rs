use crate::fem::{FeNode, FEData, FiniteElement};
use crate::fem::{FEType};
use crate::{ElementsNumbers, ElementsValues};

use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::rc::Rc;
use std::cell::RefCell;
use std::collections::HashSet;
use std::iter::FromIterator;


pub struct FEModel<T, V>
{
    pub nodes: Vec<Rc<RefCell<FeNode<T, V>>>>,
    pub elements: Vec<FiniteElement<T, V>>,
}


impl<T, V> FEModel<T, V>
    where T: Copy + PartialEq + Into<ElementsNumbers> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + From<ElementsNumbers> + Eq + Hash + SubAssign + Debug +
             Mul<Output = T> + PartialOrd + Default + Add<Output = T> + 'static,
          V: Copy + From<ElementsValues> + Sub<Output = V> + Default + Mul<Output = V> +
             Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign +
             SubAssign + Into<ElementsValues> + 'static
{
    pub fn create() -> Self
    {
        FEModel { nodes: Vec::new(), elements: Vec::new() }
    }


    pub fn add_node(&mut self, number: T, x: V, y: V, z: V) -> Result<(), String>
    {
        if self.nodes.iter().position(|node|
            node.as_ref().borrow().number == number).is_none()
        {
            let node = FeNode::create(number, x, y, z);
            self.nodes.push(Rc::new(RefCell::new(node)));
            Ok(())
        }
        else
        {
            Err(format!("FEModel: Node {} could not be added because it already exists!",
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
            self.nodes.remove(position);
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
             number already exists!", data.number.into()));
        }
        let nodes_numbers_set = HashSet::<T>::from_iter(
            nodes_numbers.iter().cloned());
        if nodes_numbers.len() != nodes_numbers_set.len()
        {
            return Err(format!("FEModel: Element {} could not be added! All nodes numbers \
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
}
