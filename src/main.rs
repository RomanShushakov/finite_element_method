use std::hash::Hash;
use std::ops::{Rem, Sub, Div, Mul, Add, SubAssign, AddAssign, MulAssign};
use std::fmt::Debug;

use extended_matrix::extended_matrix::ExtendedMatrix;

use crate::my_float::MyFloatTrait;
use crate::fem::fe_model::FEModel;
use crate::fem::finite_elements::finite_element::FEType;
use crate::fem::global_analysis::fe_boundary_condition::BCType;
use crate::fem::global_analysis::fe_dof_parameter_data::GlobalDOFParameter;

mod fem;

mod auxiliary;

mod my_float;


pub const TOLERANCE: f32 = 1e-6;


fn main() -> Result<(), String>
{
    // let mut fe_model = FEModel::<ElementsNumbers,ElementsValues>::create();
    // fe_model.add_node(1, 0.0, 4.0, 0.0)?;
    // fe_model.add_node(2, 5.0, 4.0, 0.0)?;
    // fe_model.add_node(3, 0.0, 0.0, 0.0)?;
    // fe_model.add_node(4, 3.0, 0.0, 0.0)?;
    // fe_model.add_node(5, 6.0, 0.0, 0.0)?;
    //
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![1, 2],
    //     FEData { number: 1, nodes: Vec::new(), properties: vec![2e11, 1e-6] })?;
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![1, 3],
    //     FEData { number: 2, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![2, 3],
    //     FEData { number: 3, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![2, 4],
    //     FEData { number: 4, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![2, 5],
    //     FEData { number: 5, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![3, 4],
    //     FEData { number: 6, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![4, 5],
    //     FEData { number: 7, nodes: Vec::new(), properties: vec![2e11, 1e-5, 1e-5] })?;
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 1, 1,
    //     GlobalDOFParameter::X, 0.0)?;
    // fe_model.add_bc(
    //     BCType::Displacement, 2, 1,
    //     GlobalDOFParameter::Y, 0.0)?;
    // fe_model.add_bc(
    //     BCType::Displacement, 3, 3,
    //     GlobalDOFParameter::X, 0.0)?;
    // fe_model.add_bc(
    //     BCType::Displacement, 4, 5,
    //     GlobalDOFParameter::X, 0.0)?;
    // fe_model.add_bc(
    //     BCType::Displacement, 5, 5,
    //     GlobalDOFParameter::Y, 0.0)?;
    //
    // fe_model.add_bc(
    //     BCType::Force, 1, 2,
    //     GlobalDOFParameter::X, 10000.0)?;
    // fe_model.add_bc(
    //     BCType::Force, 2, 4,
    //     GlobalDOFParameter::Y, -10000.0)?;
    //
    // let global_analysis_result = fe_model.global_analysis()?;
    // let reactions = global_analysis_result.extract_reactions();
    // for (reaction, dof_parameter_data) in
    //     reactions.reactions_values.iter().zip(reactions.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", reaction, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    // let displacements = global_analysis_result.extract_displacements();
    // for (displacement, dof_parameter_data) in
    //     displacements.displacements_values.iter().zip(displacements.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", displacement, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    //
    // fe_model.update_node(2, 3.0, 4.0, 0.0)?;
    //
    // let global_analysis_result = fe_model.global_analysis()?;
    // let reactions = global_analysis_result.extract_reactions();
    // for (reaction, dof_parameter_data) in
    //     reactions.reactions_values.iter().zip(reactions.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", reaction, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    // let displacements = global_analysis_result.extract_displacements();
    // for (displacement, dof_parameter_data) in
    //     displacements.displacements_values.iter().zip(displacements.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", displacement, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    //
    // fe_model.update_element(
    //     vec![1, 2],
    //     FEData { number: 1, nodes: Vec::new(), properties: vec![2e11, 1e-5] })?;
    //
    // let global_analysis_result = fe_model.global_analysis()?;
    // let reactions = global_analysis_result.extract_reactions();
    // for (reaction, dof_parameter_data) in
    //     reactions.reactions_values.iter().zip(reactions.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", reaction, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    // let displacements = global_analysis_result.extract_displacements();
    // for (displacement, dof_parameter_data) in
    //     displacements.displacements_values.iter().zip(displacements.dof_parameters_data.iter())
    // {
    //     println!("{}, node: {}, parameter: {:?}", displacement, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    //
    // println!();
    // let elements_analysis_results =
    //     fe_model.elements_analysis(&displacements)?;
    // println!("{:?}", elements_analysis_results);

    // let mut fe_model = FEModel::<ElementsNumbers,ElementsValues>::create();
    // fe_model.add_node(1, 0.0, 0.0, 0.0)?;
    // fe_model.add_node(2, 80.0, 0.0, 0.0)?;
    //
    //
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![1, 2],
    //     FEData { number: 1, nodes: Vec::new(), properties: vec![1.0, 1.0] })?;
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 1, 1,
    //     GlobalDOFParameter::X, 0.0)?;
    //
    // fe_model.add_bc(
    //     BCType::Force, 1, 2,
    //     GlobalDOFParameter::X, 100.0)?;
    //
    // let global_analysis_result = fe_model.global_analysis()?;
    // let reactions = global_analysis_result.extract_reactions();
    // for (reaction, dof_parameter_data) in
    //     reactions.reactions_values.iter().zip(reactions.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", reaction, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    // let displacements = global_analysis_result.extract_displacements();
    // for (displacement, dof_parameter_data) in
    //     displacements.displacements_values.iter().zip(displacements.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", displacement, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    //
    // let global_analysis_result = fe_model.global_analysis()?;
    // let reactions = global_analysis_result.extract_reactions();
    // for (reaction, dof_parameter_data) in
    //     reactions.reactions_values.iter().zip(reactions.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", reaction, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    // let displacements = global_analysis_result.extract_displacements();
    // for (displacement, dof_parameter_data) in
    //     displacements.displacements_values.iter().zip(displacements.dof_parameters_data.iter())
    // {
    //     println!("{}, node: {}, parameter: {:?}", displacement, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    //
    // println!();
    // let element_analysis_data =
    //     fe_model.elements[0].extract_element_analysis_data(&displacements)?;
    // println!("{:?}", element_analysis_data);


    let mut fe_model = FEModel::create(TOLERANCE);
    fe_model.add_node(1u32, 0f32, 0.0, 0.0)?;
    fe_model.add_node(2, 0.0, 30.0, 0.0)?;
    fe_model.add_node(3, 36.0, 0.0, 0.0)?;
    fe_model.add_node(4, 36.0, 30.0, 0.0)?;
    fe_model.add_node(5, 72.0, 0.0, 0.0)?;
    fe_model.add_node(6, 72.0, 30.0, 0.0)?;
    fe_model.add_node(7, 108.0, 0.0, 0.0)?;
    fe_model.add_node(8, 108.0, 30.0, 0.0)?;
    fe_model.add_node(9, 148.0, 0.0, 0.0)?;
    fe_model.add_node(10, 148.0, 30.0, 0.0)?;

    fe_model.add_element(1, FEType::Truss2n1ip, vec![1, 2],
        vec![1e6, 2.0])?;

    fe_model.add_element(2, FEType::Truss2n1ip, vec![1, 3],
        vec![1e6, 2.0])?;

    fe_model.add_element(3, FEType::Truss2n1ip, vec![1, 4],
        vec![1e6, 2.0])?;

    fe_model.add_element(4, FEType::Truss2n1ip, vec![2, 4],
        vec![1e6, 2.0])?;

    fe_model.add_element(5, FEType::Truss2n1ip, vec![3, 4],
        vec![1e6, 2.0])?;

    fe_model.add_element(6, FEType::Truss2n2ip, vec![3, 5],
        vec![1e6, 2.0])?;

    fe_model.add_element(7, FEType::Truss2n2ip, vec![3, 6],
        vec![1e6, 2.0])?;

    fe_model.add_element(8, FEType::Truss2n2ip, vec![4, 6],
        vec![1e6, 2.0])?;

    fe_model.add_element(9, FEType::Truss2n2ip, vec![5, 6],
        vec![1e6, 2.0])?;

    fe_model.add_element(10, FEType::Truss2n2ip, vec![5, 7],
        vec![1e6, 2.0])?;

    fe_model.add_element(11, FEType::Truss2n2ip, vec![5, 8],
        vec![1e6, 2.0])?;

    fe_model.add_element(12, FEType::Truss2n2ip, vec![6, 8],
        vec![1e6, 2.0])?;

    fe_model.add_element(13, FEType::Truss2n2ip, vec![7, 8],
        vec![1e6, 2.0])?;

    fe_model.add_element(14, FEType::Truss2n2ip, vec![7, 9],
        vec![1e6, 2.0])?;

    fe_model.add_element(15, FEType::Truss2n2ip, vec![7, 10],
        vec![1e6, 2.0])?;

    fe_model.add_element(16, FEType::Truss2n2ip, vec![8, 10],
        vec![1e6, 2.0])?;

    fe_model.add_element(17, FEType::Truss2n2ip, vec![9, 10],
        vec![1e6, 2.0])?;


    fe_model.add_bc(
        BCType::Displacement, 1, 1,
        GlobalDOFParameter::X, 0.0)?;
    fe_model.add_bc(
        BCType::Displacement, 2, 1,
        GlobalDOFParameter::Y, 0.0)?;
    fe_model.add_bc(
        BCType::Displacement, 3, 2,
        GlobalDOFParameter::X, 0.0)?;

    fe_model.add_bc(
        BCType::Force, 1, 4,
        GlobalDOFParameter::Y, 150.0)?;
    fe_model.add_bc(
        BCType::Force, 2, 6,
        GlobalDOFParameter::Y, 150.0)?;
    fe_model.add_bc(
        BCType::Force, 3, 8,
        GlobalDOFParameter::Y, 150.0)?;
    fe_model.add_bc(
        BCType::Force, 4, 10,
        GlobalDOFParameter::Y, 100.0)?;

    let global_analysis_result = fe_model.global_analysis()?;
    let reactions = global_analysis_result.extract_reactions();
    for (reaction, dof_parameter_data) in
        reactions.reactions_values.iter().zip(reactions.dof_parameters_data)
    {
        println!("{}, node: {}, parameter: {:?}", reaction, dof_parameter_data.node_number,
                 dof_parameter_data.dof_parameter);
    }
    println!();
    let displacements = global_analysis_result.extract_displacements();
    for (displacement, dof_parameter_data) in
        displacements.displacements_values.iter().zip(displacements.dof_parameters_data.iter())
    {
        println!("{}, node: {}, parameter: {:?}", displacement, dof_parameter_data.node_number,
                 dof_parameter_data.dof_parameter);
    }

    println!();
    let elements_analysis_results =
        fe_model.elements_analysis(&displacements)?;
    println!("{:?}", elements_analysis_results);



    // let mut fe_model = FEModel::create(TOLERANCE);
    // fe_model.add_node(1u16, 0.0, 0.0, 0.0)?;
    // fe_model.add_node(2, 4.0, 0.0, 0.0)?;
    // fe_model.add_node(3, 0.0, 3.0, 0.0)?;
    // fe_model.add_node(4, 4.0, 3.0, 0.0)?;
    //
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![1, 4],
    //     FEData { number: 1, nodes: Vec::new(), properties: vec![128000000.0, 0.0625] })?;
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![3, 4],
    //     FEData { number: 2, nodes: Vec::new(), properties: vec![128000000.0, 0.0625] })?;
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![2, 4],
    //     FEData { number: 3, nodes: Vec::new(), properties: vec![128000000.0, 0.0625] })?;
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 1, 1,
    //     GlobalDOFParameter::X, 0.0)?;
    // fe_model.add_bc(
    //     BCType::Displacement, 2, 1,
    //     GlobalDOFParameter::Y, 0.0)?;
    // fe_model.add_bc(
    //     BCType::Displacement, 3, 3,
    //     GlobalDOFParameter::X, 0.0)?;
    // fe_model.add_bc(
    //     BCType::Displacement, 4, 2,
    //     GlobalDOFParameter::Y, -0.025)?;
    //
    //
    //
    // let global_analysis_result = fe_model.global_analysis()?;
    // let reactions = global_analysis_result.extract_reactions();
    // for (reaction, dof_parameter_data) in
    //     reactions.reactions_values.iter().zip(reactions.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", reaction, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    // let displacements = global_analysis_result.extract_displacements();
    // for (displacement, dof_parameter_data) in
    //     displacements.displacements_values.iter().zip(displacements.dof_parameters_data.iter())
    // {
    //     println!("{}, node: {}, parameter: {:?}", displacement, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    //
    // println!();
    // let elements_analysis_results =
    //     fe_model.elements_analysis(&displacements)?;
    // println!("{:?}", elements_analysis_results);


    // let mut fe_model = FEModel::create(TOLERANCE);
    // fe_model.add_node(1u16, 0.0, 0.0, 0.0)?;
    // fe_model.add_node(2, 80.0, 0.0, 0.0)?;
    //
    // fe_model.add_element(
    //     FEType::Truss2n2ip,
    //     vec![1, 2],
    //     FEData { number: 1, nodes: Vec::new(), properties: vec![1e6, 1.0, 9.0] })?;
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 1, 1,
    //     GlobalDOFParameter::X, 0.0)?;
    //
    // fe_model.add_bc(
    //     BCType::Force, 1, 2,
    //     GlobalDOFParameter::X, 100.0)?;
    //
    // let global_analysis_result = fe_model.global_analysis()?;
    // let reactions = global_analysis_result.extract_reactions();
    // for (reaction, dof_parameter_data) in
    //     reactions.reactions_values.iter().zip(reactions.dof_parameters_data)
    // {
    //     println!("{}, node: {}, parameter: {:?}", reaction, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    // println!();
    // let displacements = global_analysis_result.extract_displacements();
    // for (displacement, dof_parameter_data) in
    //     displacements.displacements_values.iter().zip(displacements.dof_parameters_data.iter())
    // {
    //     println!("{}, node: {}, parameter: {:?}", displacement, dof_parameter_data.node_number,
    //              dof_parameter_data.dof_parameter);
    // }
    //
    // println!();
    // let elements_analysis_results =
    //     fe_model.elements_analysis(&displacements)?;
    // println!("{:?}", elements_analysis_results);

    Ok(())
}
