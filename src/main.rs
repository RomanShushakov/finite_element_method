mod extended_matrix;
use extended_matrix::ExtendedMatrix;
use extended_matrix::{return_symmetric_matrix_struct, return_non_symmetric_matrix_struct};

mod fem;
use crate::fem::{FeNode, GlobalCoordinates, FEModel, FEType, FEData};
use fem::finite_elements::truss::truss2n2ip::Truss2n2ip;
use crate::extended_matrix::basic_matrix::basic_matrix::BasicMatrixTrait;

use std::mem;


pub type ElementsNumbers = u16;
pub type ElementsValues = f32;


pub const TOLERANCE: ElementsValues = 1e-6;


fn main() -> Result<(), String>
{
    let mut fe_model = FEModel::<ElementsNumbers,_>::create();
    fe_model.add_node(1, 4.0, 0.0, 0.0)?;
    fe_model.add_node(2, 4.0, 3.0, 0.0)?;
    fe_model.add_node(3, 0.0, 0.0, 0.0)?;
    fe_model.add_node(4, 0.0, 3.0, 0.0)?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![1, 2],
        FEData { number: 1, nodes: Vec::new(), properties: vec![128000000.0, 0.0625] })?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![2, 3],
        FEData { number: 2, nodes: Vec::new(), properties: vec![128000000.0, 0.0625] })?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![2, 4],
        FEData { number: 3, nodes: Vec::new(), properties: vec![128000000.0, 0.0625] })?;
    let mut global_stiffness_matrix = ExtendedMatrix::create(
        24 as ElementsNumbers,
        24 as ElementsNumbers,
        vec![0.0; 24 * 24]);
    for element in &fe_model.elements
    {
        let element_stiffness_matrix = element.extract_stiffness_matrix()?;
        let element_stiffness_groups = element.extract_stiffness_groups();
        for element_stiffness_group in element_stiffness_groups
        {
            if let Some(position) = fe_model.stiffness_groups
                .iter()
                .position(|group|
                              { group.stiffness_type == element_stiffness_group.stiffness_type &&
                                group.number_1 == element_stiffness_group.number_1 &&
                                group.number_2 == element_stiffness_group.number_2 })
            {

                global_stiffness_matrix.add_sub_matrix(
                    &element_stiffness_matrix,
                    fe_model.stiffness_groups[position].positions.as_slice(),
                    element_stiffness_group.positions.as_slice());

            }
        }


    }
    global_stiffness_matrix.remove_zeros_rows_columns();
    global_stiffness_matrix.show_matrix();


    // let m = fe_model.elements[0].extract_stiffness_matrix()?;
    // m.show_matrix();
    // println!("{:?}", m.basic_matrix.define_type());
    // println!();
    // println!("{:?}", fe_model.nodes);
    // println!();
    //
    // fe_model.update_node(2, 0.0, 3.0, 0.0)?;
    // let m = fe_model.elements[0].extract_stiffness_matrix()?;
    // m.show_matrix();
    // println!("{:?}", m.basic_matrix.define_type());
    // println!();
    // println!("{:?}", fe_model.elements[0].extract_stiffness_groups());
    // println!();
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![3, 1],
    //     FEData { number: 2, nodes: Vec::new(), properties: vec![1e6, 1.0, 9.0] })?;
    // let m = fe_model.elements[1].extract_stiffness_matrix()?;
    // m.show_matrix();
    // println!();
    // fe_model.delete_element(1)?;
    // fe_model.update_element(
    //     vec![3, 2],
    //     FEData { number: 2, nodes: Vec::new(), properties: vec![1.6e6, 3.0, 9.0] })?;
    // let m = fe_model.elements[0].extract_stiffness_matrix()?;
    // m.show_matrix();
    // println!();
    // println!("{}, {}", fe_model.nodes.len(), fe_model.elements.len());
    // println!();
    // println!("{:?}", fe_model.stiffness_groups);
    // println!();
    // fe_model.delete_node(1)?;
    // println!("{:?}", fe_model.stiffness_groups);
    // println!();
    // fe_model.delete_node(3)?;
    // println!("{:?}", fe_model.stiffness_groups);
    // println!();
    Ok(())
}
