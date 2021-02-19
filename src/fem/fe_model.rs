use crate::fem::{FeNode, Truss2n2ip};
use crate::{ElementsNumbers, ElementsValues};

use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;
use crate::extended_matrix::One;


pub struct FEModel<'a, T, V>
{
    pub nodes: Vec<FeNode<T, V>>,
    pub elements: Vec<Truss2n2ip<'a, T, V>>,
}


impl<'a, T, V> FEModel<'a, T, V>
    where T: Copy + PartialEq + Into<ElementsNumbers> + Sub<Output = T> + Div<Output = T> +
             Rem<Output = T> + From<ElementsNumbers> + Eq + Hash + SubAssign + Debug +
             Mul<Output = T> + PartialOrd + Default + Add<Output = T> + 'static,
          V: Copy + From<ElementsValues> + Sub<Output = V> + Default + Mul<Output = V> +
             Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign + MulAssign +
             SubAssign + One + Into<ElementsValues> + 'static
{
    pub fn create() -> Self
    {
        FEModel { nodes: Vec::new(), elements: Vec::new() }
    }


    pub fn add_node(&mut self, node: FeNode<T, V>)
    {
        self.nodes.push(node);
    }


    pub fn update_node(&mut self, number: T, x: V, y: V, z: V) -> Result<(), String>
    {
        if let Some(position) = self.nodes.iter().position(|node|
            node.number == number)
        {
            self.nodes[position].update(x, y, z);
            return Ok(());
        }
        Err(format!("FEModel: Node {} could not be updated!", number.into()))
    }


    pub fn add_element(&'a mut self, element_number: T, node_1_number: T, node_2_number: T,
        young_modulus: V, area: V, area_2: Option<V>) -> Result<(), String>
    {
        if let Some(pos_1) = self.nodes.iter().position(|node|
            node.number == node_1_number)
        {
            if let Some(pos_2) = self.nodes.iter().position(|node|
                node.number == node_2_number)
            {

                let element = Truss2n2ip::create(
                    element_number, &self.nodes[pos_1], &self.nodes[pos_2],
                    young_modulus, area, area_2)?;
                self.elements.push(element);
            }
            else
            {
                return Err(format!("FEModel: Element {} could not be added! Node {} does not exist!",
                               element_number.into(), node_2_number.into()));
            }
        }
        else
        {
            return Err(format!("FEModel: Element {} could not be added! Node {} does not exist!",
                               element_number.into(), node_1_number.into()));
        }
        Ok(())
    }
}
