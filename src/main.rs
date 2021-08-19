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
    // let mut fe_model = FEModel::create(TOLERANCE);
    // fe_model.add_node(1u32, 0f32, 0.0, 0.0)?;
    // fe_model.add_node(2, 0.0, 30.0, 0.0)?;
    // fe_model.add_node(3, 36.0, 0.0, 0.0)?;
    // fe_model.add_node(4, 36.0, 30.0, 0.0)?;
    // fe_model.add_node(5, 72.0, 0.0, 0.0)?;
    // fe_model.add_node(6, 72.0, 30.0, 0.0)?;
    // fe_model.add_node(7, 108.0, 0.0, 0.0)?;
    // fe_model.add_node(8, 108.0, 30.0, 0.0)?;
    // fe_model.add_node(9, 148.0, 0.0, 0.0)?;
    // fe_model.add_node(10, 148.0, 30.0, 0.0)?;
    //
    // fe_model.add_element(1, FEType::Truss2n1ip, vec![1, 2],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(2, FEType::Truss2n1ip, vec![1, 3],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(3, FEType::Truss2n1ip, vec![1, 4],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(4, FEType::Truss2n1ip, vec![2, 4],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(5, FEType::Truss2n1ip, vec![3, 4],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(6, FEType::Truss2n2ip, vec![3, 5],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(7, FEType::Truss2n2ip, vec![3, 6],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(8, FEType::Truss2n2ip, vec![4, 6],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(9, FEType::Truss2n2ip, vec![5, 6],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(10, FEType::Truss2n2ip, vec![5, 7],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(11, FEType::Truss2n2ip, vec![5, 8],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(12, FEType::Truss2n2ip, vec![6, 8],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(13, FEType::Truss2n2ip, vec![7, 8],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(14, FEType::Truss2n2ip, vec![7, 9],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(15, FEType::Truss2n2ip, vec![7, 10],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(16, FEType::Truss2n2ip, vec![8, 10],
    //     vec![1e6, 2.0])?;
    //
    // fe_model.add_element(17, FEType::Truss2n2ip, vec![9, 10],
    //     vec![1e6, 2.0])?;
    //
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 1, 1,
    //     GlobalDOFParameter::X, 0.0)?;
    // fe_model.add_bc(
    //     BCType::Displacement, 2, 1,
    //     GlobalDOFParameter::Y, 0.0)?;
    // fe_model.add_bc(
    //     BCType::Displacement, 3, 2,
    //     GlobalDOFParameter::X, 0.0)?;
    //
    // fe_model.add_bc(
    //     BCType::Force, 1, 4,
    //     GlobalDOFParameter::Y, 150.0)?;
    // fe_model.add_bc(
    //     BCType::Force, 2, 6,
    //     GlobalDOFParameter::Y, 150.0)?;
    // fe_model.add_bc(
    //     BCType::Force, 3, 8,
    //     GlobalDOFParameter::Y, 150.0)?;
    // fe_model.add_bc(
    //     BCType::Force, 4, 10,
    //     GlobalDOFParameter::Y, 100.0)?;
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
    // fe_model.add_node(2, 4.0, 0.0, 0.0)?;
    // fe_model.add_node(3, 0.0, 3.0, 0.0)?;
    // fe_model.add_node(4, 4.0, 3.0, 0.0)?;
    //
    // fe_model.add_element(
    //     1,
    //     FEType::Truss2n1ip,
    //     vec![1, 4],
    //     vec![128000000.0, 0.0625])?;
    // fe_model.add_element(
    //     2,
    //     FEType::Truss2n1ip,
    //     vec![3, 4],
    //     vec![128000000.0, 0.0625])?;
    // fe_model.add_element(
    //     3,
    //     FEType::Truss2n1ip,
    //     vec![2, 4],
    //     vec![128000000.0, 0.0625])?;
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
    // let global_analysis_result = fe_model.global_analysis()?;
    // let reactions = global_analysis_result.extract_reactions();
    // for (reaction, dof_parameter_data) in
    //     reactions.extract_reactions_values().iter()
    //         .zip(reactions.extract_dof_parameters_data())
    // {
    //     println!("{}, node: {}, parameter: {:?}", reaction, dof_parameter_data.extract_node_number(),
    //              dof_parameter_data.extract_dof_parameter());
    // }
    // println!();
    // let displacements = global_analysis_result.extract_displacements();
    // for (displacement, dof_parameter_data) in
    //     displacements.extract_displacements_values().iter()
    //         .zip(displacements.extract_dof_parameters_data().iter())
    // {
    //     println!("{}, node: {}, parameter: {:?}", displacement, dof_parameter_data.extract_node_number(),
    //              dof_parameter_data.extract_dof_parameter());
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
    //     1,
    //     FEType::Truss2n1ip,
    //     vec![1, 2],
    //     vec![1e6, 1.0, 9.0])?;
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



    let mut fe_model = FEModel::create(TOLERANCE);
    fe_model.add_node(1u16, 0.0, 0.0, 0.0)?;
    fe_model.add_node(2, 25.0, 0.0, 0.0)?;
    fe_model.add_node(3, 50.0, 0.0, 0.0)?;
    fe_model.add_node(4, 75.0, 0.0, 0.0)?;
    fe_model.add_node(5, 100.0, 0.0, 0.0)?;

    fe_model.add_element(
        1,
        FEType::Beam2n1ipT,
        vec![1, 2],
        vec![1e6, 0.3, 1.0, 2.0, 1.0, 0.0, 1.0, 0.8333, 0.0, 0.0, -1.0])?;

    fe_model.add_element(
        2,
        FEType::Beam2n1ipT,
        vec![2, 3],
        vec![1e6, 0.3, 1.0, 2.0, 1.0, 0.0, 1.0, 0.8333, 0.0, 0.0, -1.0])?;

    fe_model.add_element(
        3,
        FEType::Beam2n1ipT,
        vec![3, 4],
        vec![1e6, 0.3, 1.0, 2.0, 1.0, 0.0, 1.0, 0.8333, 0.0, 0.0, -1.0])?;

    fe_model.add_element(
        4,
        FEType::Beam2n1ipT,
        vec![4, 5],
        vec![1e6, 0.3, 1.0, 2.0, 1.0, 0.0, 1.0, 0.8333, 0.0, 0.0, -1.0])?;

    fe_model.add_bc(
        BCType::Displacement, 1, 1,
        GlobalDOFParameter::X, 0.0)?;

    fe_model.add_bc(
        BCType::Displacement, 2, 1,
        GlobalDOFParameter::Y, 0.0)?;

    fe_model.add_bc(
        BCType::Displacement, 3, 1,
        GlobalDOFParameter::Z, 0.0)?;

    fe_model.add_bc(
        BCType::Displacement, 4, 1,
        GlobalDOFParameter::ThX, 0.0)?;

    fe_model.add_bc(
        BCType::Displacement, 5, 1,
        GlobalDOFParameter::ThY, 0.0)?;

    fe_model.add_bc(
        BCType::Displacement, 6, 1,
        GlobalDOFParameter::ThZ, 0.0)?;

    // fe_model.add_bc(
    //     BCType::Displacement, 7, 2,
    //     GlobalDOFParameter::X, 0.0)?;
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 8, 2,
    //     GlobalDOFParameter::Y, 0.0)?;
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 9, 2,
    //     GlobalDOFParameter::Z, 0.0)?;
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 10, 2,
    //     GlobalDOFParameter::ThX, 0.0)?;
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 11, 2,
    //     GlobalDOFParameter::ThY, 0.0)?;
    //
    // fe_model.add_bc(
    //     BCType::Displacement, 12, 2,
    //     GlobalDOFParameter::ThZ, 0.0)?;

    fe_model.add_bc(
        BCType::Force, 1, 5,
        GlobalDOFParameter::Y, 100.0)?;

    let global_analysis_result = fe_model.global_analysis()?;
    let reactions = global_analysis_result.extract_reactions();
    for (reaction, dof_parameter_data) in
        reactions.extract_reactions_values().iter()
            .zip(reactions.extract_dof_parameters_data())
    {
        println!("reaction: {}, node: {}, parameter: {:?}", reaction, dof_parameter_data.extract_node_number(),
                 dof_parameter_data.extract_dof_parameter());
    }
    println!();
    let displacements = global_analysis_result.extract_displacements();
    for (displacement, dof_parameter_data) in
        displacements.extract_displacements_values().iter()
            .zip(displacements.extract_dof_parameters_data().iter())
    {
        println!("displacement: {}, node: {}, parameter: {:?}", displacement, dof_parameter_data.extract_node_number(),
                 dof_parameter_data.extract_dof_parameter());
    }

    println!();
    let elements_analysis_results =
        fe_model.elements_analysis(&displacements)?;
    println!("{:?}", elements_analysis_results);

    Ok(())
}
