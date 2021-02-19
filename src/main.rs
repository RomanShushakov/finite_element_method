mod extended_matrix;
use extended_matrix::ExtendedMatrix;
use extended_matrix::{return_symmetric_matrix_struct, return_non_symmetric_matrix_struct};

mod fem;
use crate::fem::{FeNode, GlobalCoordinates, FEModel};
use fem::finite_elements::truss::truss2n2ip::Truss2n2ip;
use crate::extended_matrix::basic_matrix::basic_matrix::BasicMatrix;

use std::mem;


pub type ElementsNumbers = u16;
pub type ElementsValues = f64;


pub const TOLERANCE: ElementsValues = 1e-9;


fn main() -> Result<(), String>
{


    let mut fe_model = FEModel::create();
    let node_1 = FeNode::create(1u16, 0.0, 0.0, 0.0);
    fe_model.add_node(node_1);
    let node_2 = FeNode::create(2u16, 4.0, 3.0, 0.0);
    fe_model.add_node(node_2);
    let node_3 = FeNode::create(3u16, 80.0, 0.0, 0.0);
    fe_model.add_node(node_3);

    fe_model.add_element(1u16, 1u16, 2u16,
    128000000.0, 0.0625, None)?;

    // fe_model.add_element(2u16, 10u16, 2u16,
    // 128000000.0, 0.0625, None)?;


    // let mut elem_1 = Truss2n2ip::create(
    // 1u16, &fe_model.nodes[0], &fe_model.nodes[1],
    // 128000000.0, 0.0625, Some(0.0625)).unwrap();
    // fe_model.elements = Some(vec![elem_1]);

    // fe_model.update_node(3u16, 0.0, 3.0, 0.0)?;
    //
    let node_4 = FeNode::create(4u16, 3.0, 3.0, 3.0);
    fe_model.add_node(node_4);
    //
    // println!("{:?}", fe_model.nodes);


    // match fe_model.update_node(2u16, 0.0, 3.0, 0.0)
    // {
    //     Ok(_) => println!("Ok"),
    //     Err(e) => println!("{}", e),
    // }

    // if let Some(mut elements) = fe_model.elements
    // {
    //     elements[0].state.rotation_matrix.show_matrix();
    //     println!();
    //     // elements[0].update(&fe_model.nodes[0], &fe_model.nodes[2],
    //     // 128000000.0, 0.0625, Some(0.0625));
    //     // elements[0].state.rotation_matrix.show_matrix();
    // }
    //


    Ok(())

    // fe_model.elements.unwrap()[0].update(&fe_model.nodes[0], &fe_model.nodes[2],
    // 128000000.0, 0.0625, Some(0.0625));
    // elem.state.local_stiffness_matrix.show_matrix();
    // println!();
    // let matrix = elem.extract_stiffness_matrix().unwrap();
    // matrix.show_matrix();
    // println!();
    // println!("{:?}", matrix.basic_matrix.define_type());
    // println!();
    // let stiffness_groups = elem.extract_stiffness_groups();
    // println!("{:?}", stiffness_groups);

}
