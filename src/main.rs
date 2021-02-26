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
    fe_model.add_node(1, 0.0, 4.0, 0.0)?;
    fe_model.add_node(2, 3.0, 4.0, 0.0)?;
    fe_model.add_node(3, 0.0, 0.0, 0.0)?;
    fe_model.add_node(4, 3.0, 0.0, 0.0)?;
    fe_model.add_node(5, 6.0, 0.0, 0.0)?;

    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![1, 2],
        FEData { number: 1, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![1, 3],
        FEData { number: 2, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![2, 3],
        FEData { number: 3, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![2, 4],
        FEData { number: 4, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![5, 2],
        FEData { number: 5, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![3, 4],
        FEData { number: 6, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    fe_model.add_element(
        FEType::Truss2n2ip,
        vec![4, 5],
        FEData { number: 7, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;

    fe_model.add_bc(
        BCType::Displacement, 1, 1,
        GlobalDOFParameter::X, 0.0)?;
    fe_model.add_bc(
        BCType::Displacement, 2, 1,
        GlobalDOFParameter::Y, 0.0)?;
    fe_model.add_bc(
        BCType::Displacement, 3, 3,
        GlobalDOFParameter::X, 0.0)?;
    fe_model.add_bc(
        BCType::Displacement, 4, 5,
        GlobalDOFParameter::X, 0.0)?;
    fe_model.add_bc(
        BCType::Displacement, 5, 5,
        GlobalDOFParameter::Y, 0.0)?;

    fe_model.add_bc(
        BCType::Force, 1, 2,
        GlobalDOFParameter::X, 10000.0)?;
    fe_model.add_bc(
        BCType::Force, 2, 4,
        GlobalDOFParameter::Y, -10000.0)?;

    let global_analysis_result = fe_model.global_analysis()?;
    global_analysis_result.show_reactions();
    println!();
    global_analysis_result.show_displacements();

    Ok(())
}
