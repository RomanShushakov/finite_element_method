// use std::fmt::Debug;

use extended_matrix::FloatTrait;

use crate::fem::FEM;
use crate::fem::structs::Truss;


enum TrussError
{
    Number(u32),
    SameNodes(u32, u32),
}


impl TrussError
{
    fn compose_error_message(&self) -> String
    {
        match self
        {
            TrussError::Number(number) => format!("Truss element with number {number} already exists!"),
            TrussError::SameNodes(node_1_number, node_2_number) => 
            {
                format!("Truss element with node number {node_1_number} and {node_2_number} already exists!")
            },
        }
    }
}


impl<V> FEM<V>
    where V: FloatTrait<Output = V>
{
    fn check_truss_data(&self, number: u32, node_1_number: u32, node_2_number: u32) -> Option<TrussError>
    {
        for (truss_number, truss_element) in self.get_truss_elements().iter()
        {
            if *truss_number == number
            {
                return Some(TrussError::Number(number));
            }
            if truss_element.is_nodes_numbers_same(node_1_number, node_2_number)
            {
                return Some(TrussError::SameNodes(node_1_number, node_2_number));
            }
        }
        None
    }


    // pub fn add_node(&mut self, number: u32, x: V, y: V, z: V) -> Result<(), String>
    // {
    //     if let Some(node_error) = self.check_node_data(number, *self.get_nodes_count(), x, y, z)
    //     {
    //         return Err(node_error.compose_error_message());
    //     }

    //     let node = Node::create(*self.get_nodes_count(), x, y, z);
    //     self.get_mut_nodes().insert(number, node);

    //     *self.get_mut_nodes_count() += 1;
    //     Ok(())
    // }

    pub fn add_truss(&mut self,
        number: u32,
        node_1_number: u32,
        node_2_number: u32,
        young_modulus: V,
        poisson_ratio: V,
        area: V,
        optional_area_2: Option<V>,
    )
        -> Result<(), String>
    {
        self.check_node_exist(node_1_number)?;
        self.check_node_exist(node_2_number)?;
            if let Some(truss_error) = self.check_truss_data(number, node_1_number, node_2_number)
        {
            return Err(truss_error.compose_error_message());
        }

        let truss_element = Truss::create(
            node_1_number, 
            node_2_number, 
            young_modulus, 
            poisson_ratio, 
            area, 
            optional_area_2, 
            self.get_nodes(), 
            self.get_props().get_rel_tol(),
            self.get_props().get_abs_tol(),
        )?;
        self.get_mut_truss_elements().insert(number, truss_element);

        Ok(())
    }
}
