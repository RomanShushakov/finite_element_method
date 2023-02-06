use std::fmt::Debug;

use extended_matrix::FloatTrait;

use crate::fem::FEM;
use crate::fem::structs::Node;


enum NodeError<V>
{
    Number(u32),
    NumberNotExist(u32),
    Index(usize),
    Coordinates(V, V, V),
}


impl<V> NodeError<V>
    where V: Debug
{
    fn compose_error_message(&self) -> String
    {
        match self
        {
            NodeError::Number(number) => format!("Node with number {number} already exists!"),
            NodeError::NumberNotExist(number) => format!("Node with number {number} does not exist!"),
            NodeError::Index(index) => format!("Node wit index {index} already exists!"),
            NodeError::Coordinates(x, y, z) => 
            {
                format!("Node with coordinates x: {x:?}, y: {y:?}, z: {z:?} already exists!")
            }
        }
    }
}


impl<V> FEM<V>
    where V: FloatTrait
{
    fn check_node_data(&self, number: u32, index: usize, x: V, y: V, z: V) -> Option<NodeError<V>>
    {
        for (node_number, node) in self.get_nodes().iter()
        {
            if *node_number == number
            {
                return Some(NodeError::Number(number));
            }
            if node.is_index_same(index)
            {
                return Some(NodeError::Index(index));
            }
            if node.is_coordinates_same(x, y, z)
            {
                return Some(NodeError::Coordinates(x, y, z));
            }
        }
        None
    }


    pub fn add_node(&mut self, number: u32, x: V, y: V, z: V) -> Result<(), String>
    {
        if let Some(node_error) = self.check_node_data(number, *self.get_index(), x, y, z)
        {
            return Err(node_error.compose_error_message());
        }

        let node = Node::create(*self.get_index(), x, y, z);
        self.get_mut_nodes().insert(number, node);

        *self.get_mut_index() += 1;
        Ok(())
    }


    pub(crate) fn check_node_exist(&self, number: u32) -> Result<(), String>
    {
        if !self.get_nodes().contains_key(&number) 
        {
            return Err(NodeError::<V>::NumberNotExist(number).compose_error_message());
        }
        Ok(())
    }
}
