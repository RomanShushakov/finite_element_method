use crate::fem::{FeNode, FEData, FiniteElement, StiffnessGroup, StiffnessType};
use crate::fem::{FEType};
use crate::fem::{STIFFNESS_TYPES_NUMBER};
use crate::{ElementsNumbers, ElementsValues};

use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use std::rc::Rc;
use std::cell::RefCell;
use std::collections::HashSet;
use std::iter::FromIterator;


const GLOBAL_DOF: ElementsNumbers = 6;


pub struct FEModel<T, V>
{
    pub nodes: Vec<Rc<RefCell<FeNode<T, V>>>>,
    pub elements: Vec<FiniteElement<T, V>>,
    pub stiffness_groups: Vec<StiffnessGroup<T>>
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
        FEModel { nodes: Vec::new(), elements: Vec::new(), stiffness_groups: Vec::new() }
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
            let mut position = 0;
            let columns_number = nodes_numbers.len();
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
                        let row = position / columns_number;
                        let column = position % columns_number;
                        println!("{}, {}", row, column);
                        for k in 0..STIFFNESS_TYPES_NUMBER
                        {
                            for m in (row * 6 + k / 2 * 3)..(row * 6 + k / 2 * 3 + 3)
                            {
                                for n in (column * 6 + k % 2 * 3)..(column * 6 + k % 2 * 3 + 3)
                                {
                                    print!("{}, {}; ", m, n);
                                }
                            }
                            println!();
                            let stiffness_type = StiffnessType::iterator()
                                .nth(k)
                                .ok_or("FEModel: Stiffness type could not be defined")?;
                            let stiffness_group = StiffnessGroup { stiffness_type: *stiffness_type,
                                number_1: nodes_numbers[j],
                                number_2: nodes_numbers[j],
                                positions: Vec::new(),
                            };
                            stiffness_groups.push(stiffness_group);
                        }
                        println!();
                        position += 1;
                    }
                    let row = position / columns_number;
                    let column = position % columns_number;
                    println!("{}, {}", row, column);
                    for k in 0..STIFFNESS_TYPES_NUMBER
                    {
                        for m in (row * 6 + k / 2 * 3)..(row * 6 + k / 2 * 3 + 3)
                        {
                            for n in (column * 6 + k % 2 * 3)..(column * 6 + k % 2 * 3 + 3)
                            {
                                print!("{}, {}; ", m, n);
                            }
                        }
                        println!();
                        let stiffness_type = StiffnessType::iterator()
                            .nth(k)
                            .ok_or("FEModel: Stiffness type could not be defined")?;
                        let stiffness_group = StiffnessGroup { stiffness_type: *stiffness_type,
                            number_1: excluded,
                            number_2: v_lhs[j],
                            positions: Vec::new(),
                        };
                        stiffness_groups.push(stiffness_group);
                    }
                    println!();
                    position += 1;
                }
            }
            for i in 0..nodes_numbers.len() - 1
            {
                let row = position / columns_number;
                let column = position % columns_number;
                println!("{}, {}", row, column);
                for k in 0..STIFFNESS_TYPES_NUMBER
                {
                    for m in (row * 6 + k / 2 * 3)..(row * 6 + k / 2 * 3 + 3)
                    {
                        for n in (column * 6 + k % 2 * 3)..(column * 6 + k % 2 * 3 + 3)
                        {
                            print!("{}, {}; ", m, n);
                        }
                    }
                    println!();
                    let stiffness_type = StiffnessType::iterator()
                        .nth(k)
                        .ok_or("FEModel: Stiffness type could not be defined")?;
                    let stiffness_group = StiffnessGroup { stiffness_type: *stiffness_type,
                        number_1: nodes_numbers[nodes_numbers.len() - 1],
                        number_2: nodes_numbers[i],
                        positions: Vec::new(),
                    };
                    stiffness_groups.push(stiffness_group);
                }
                position += 1;
            }
            let row = position / columns_number;
            let column = position % columns_number;
            println!();
            println!("{}, {}", row, column);
            for k in 0..STIFFNESS_TYPES_NUMBER
            {
                for m in (row * 6 + k / 2 * 3)..(row * 6 + k / 2 * 3 + 3)
                {
                    for n in (column * 6 + k % 2 * 3)..(column * 6 + k % 2 * 3 + 3)
                    {
                        print!("{}, {}; ", m, n);
                    }
                }
                println!();
                let stiffness_type = StiffnessType::iterator()
                    .nth(k)
                    .ok_or("FEModel: Stiffness type could not be defined")?;
                let stiffness_group = StiffnessGroup { stiffness_type: *stiffness_type,
                    number_1: nodes_numbers[nodes_numbers.len() - 1],
                    number_2: nodes_numbers[nodes_numbers.len() - 1],
                    positions: Vec::new(),
                };
                stiffness_groups.push(stiffness_group);
            }
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
