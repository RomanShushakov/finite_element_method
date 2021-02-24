mod extended_matrix;
use extended_matrix::ExtendedMatrix;
use extended_matrix::{return_symmetric_matrix_struct, return_non_symmetric_matrix_struct};

mod fem;
use crate::fem::
    {
        FeNode, GlobalCoordinates, FEModel, FEType, FEData, GLOBAL_DOF, DOFParameterData,
        GlobalDOFParameter
    };
use fem::finite_elements::truss::truss2n2ip::Truss2n2ip;
use crate::extended_matrix::basic_matrix::basic_matrix::BasicMatrixTrait;

use std::mem;
use crate::extended_matrix::MatrixElementPosition;


pub type ElementsNumbers = u16;
pub type ElementsValues = f32;



pub const TOLERANCE: ElementsValues = 1e-6;


fn main() -> Result<(), String>
{
    let mut fe_model = FEModel::<ElementsNumbers,ElementsValues>::create();
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
    let mut global_stiffness_matrix =
        fe_model.compose_global_stiffness_matrix()?;
    global_stiffness_matrix.show_matrix();
    println!();
    let removed_zeros_rows_columns =
        global_stiffness_matrix.remove_zeros_rows_columns();
    for removed_zeros_row_column in &removed_zeros_rows_columns
    {
        let removed_force_node_number =
            fe_model.nodes[(removed_zeros_row_column.row / GLOBAL_DOF) as usize]
                .as_ref().borrow().number;
        let removed_force_component =
            GlobalDOFParameter::iterator()
                .nth((removed_zeros_row_column.row % GLOBAL_DOF) as usize)
                .ok_or("Could not find force component!")?;
        println!("{} {:?}", removed_force_node_number, removed_force_component);

    }
    println!();
    global_stiffness_matrix.show_matrix();
    println!();
    // println!("{:?}", removed_zeros_rows_columns);
    // println!();
    let separated_matrix = global_stiffness_matrix.separate(
        vec![
            MatrixElementPosition { row: 1, column: 1 },
            MatrixElementPosition { row: 4, column: 4 }])?;
    separated_matrix.k_aa.show_matrix();
    println!();
    separated_matrix.k_ab.show_matrix();
    println!();
    separated_matrix.k_ba.show_matrix();
    println!();
    separated_matrix.k_bb.show_matrix();
    println!();


    // fe_model.add_load(
    //     1, 1, GlobalDOFParameter::X, 1000.0)?;
    // println!("{:?}", fe_model.applied_loads);
    // println!();
    // fe_model.add_load(
    //     5, 2, GlobalDOFParameter::ThX, 1500.0)?;
    // println!("{:?}", fe_model.applied_loads);
    // println!();
    // fe_model.update_load(
    //     5, 4, GlobalDOFParameter::ThX, 1500.0)?;
    // println!("{:?}", fe_model.applied_loads);
    // println!();
    // fe_model.delete_load(1)?;
    // println!("{:?}", fe_model.applied_loads);
    // println!();
    //
    // fe_model.add_displacement(
    //     1, 1, GlobalDOFParameter::X, 1000.0)?;
    // println!("{:?}", fe_model.applied_displacements);
    // println!();
    // fe_model.add_displacement(
    //     5, 2, GlobalDOFParameter::ThX, 1500.0)?;
    // println!("{:?}", fe_model.applied_displacements);
    // println!();
    // fe_model.update_displacement(
    //     5, 4, GlobalDOFParameter::ThX, 1500.0)?;
    // println!("{:?}", fe_model.applied_displacements);
    // println!();
    // fe_model.delete_displacement(1)?;
    // println!("{:?}", fe_model.applied_displacements);
    // println!();
    //
    // fe_model.delete_node(1)?;
    // println!("{:?}", fe_model.applied_displacements);
    // println!();
    // println!("{:?}", fe_model.applied_loads);
    // println!();
    // println!("{:?}", fe_model.elements.len());
    // println!();

    Ok(())
}
