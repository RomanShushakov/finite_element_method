mod extended_matrix;
use extended_matrix::ExtendedMatrix;
use extended_matrix::{return_symmetric_matrix_struct, return_non_symmetric_matrix_struct};

mod fem;
use crate::fem::{FeNode, GlobalCoordinates, FEModel, FEType, FEData, GLOBAL_DOF, DOFParameterData, GlobalDOFParameter, BCType};
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
        vec![3, 2],
        FEData { number: 2, nodes: Vec::new(), properties: vec![128000000.0, 0.0625] })?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![2, 4],
        FEData { number: 3, nodes: Vec::new(), properties: vec![128000000.0, 0.0625] })?;

    fe_model.add_bc(
        BCType::Displacement, 1, 1,
        GlobalDOFParameter::Y, -0.025)?;
    fe_model.add_bc(
        BCType::Displacement, 2, 3,
        GlobalDOFParameter::X, 0.0)?;
    fe_model.add_bc(
        BCType::Displacement, 3, 3,
        GlobalDOFParameter::Y, 0.0)?;
    fe_model.add_bc(
        BCType::Displacement, 4, 4,
        GlobalDOFParameter::X, 0.0)?;

    fe_model.global_analysis()?;

    // let mut global_stiffness_matrix =
    //     fe_model.compose_global_stiffness_matrix()?;
    // global_stiffness_matrix.show_matrix();
    // println!();
    // let removed_zeros_rows_columns =
    //     global_stiffness_matrix.remove_zeros_rows_columns();
    // println!("{:?}", removed_zeros_rows_columns);
    // println!();
    // global_stiffness_matrix.show_matrix();
    // println!();
    // let separated_matrix = global_stiffness_matrix.separate(
    //     vec![
    //         MatrixElementPosition { row: 1, column: 1 },
    //         MatrixElementPosition { row: 4, column: 4 }])?;
    // separated_matrix.k_aa.show_matrix();
    // println!();
    // separated_matrix.k_ab.show_matrix();
    // println!();
    // separated_matrix.k_ba.show_matrix();
    // println!();
    // separated_matrix.k_bb.show_matrix();
    // println!();
    //
    // println!("{:?}", fe_model.state.nodes_dof_parameters_global);

    Ok(())
}
