mod extended_matrix;
use extended_matrix::extended_matrix::{ExtendedMatrix};
use extended_matrix::basic_matrix::return_symmetric_matrix_struct;
use extended_matrix::basic_matrix::return_non_symmetric_matrix_struct;
mod fem;
use fem::finite_elements::fe_node::{FeNode, GlobalCoordinates};
use fem::finite_elements::truss::truss2n2ip::Truss2n2ip;


use std::mem;
use crate::extended_matrix::basic_matrix::BasicMatrix;



pub type ElementsNumbers = u16;
pub type ElementsValues = f64;


pub const TOLERANCE: ElementsValues = 1e-9;


fn main()
{
    let m_10 = ExtendedMatrix::create(3u16, 3u16, vec![1.0, 0.0, 0.0, 0.03333333333333333, 1.0, 0.0, 0.09999999999999999, -0.027129938124702525, 1.0]);
    println!("{:?}\n", m_10.basic_matrix.define_type());
    let m_11 = ExtendedMatrix::create(3u16, 3u16, vec![3.0, -0.1, -0.2, 0.0, 7.003333333333333, -0.29333333333333333, 0.0, 0.0, 10.012041884816753]);
    let m_12 = &m_10.multiply_by_matrix(&m_11);
    match m_12
    {
        Ok(matrix) =>
            {
                matrix.show_matrix();
                // let basic_matrix = return_symmetric_matrix_struct(matrix.basic_matrix.clone());
                // println!("{:?}\n", basic_matrix);
                // println!("{:?}\n", mem::size_of_val(&basic_matrix));
            },
        Err(e) => println!("{}\n", e),
    }
    println!();

    let m_13 = ExtendedMatrix::create(
        3u16, 3u16,
        vec![3.0, -0.1, -0.2, 0.1, 7.0, -0.3, 0.3, -0.2, 10.0]);

    let m_14 = ExtendedMatrix::create(
        3u16, 1u16, vec![7.85, -19.3, 71.4]);
    if let Ok(m_15) = &m_13.naive_gauss_elimination(&m_14)
    {
        m_15.show_matrix();
    }
    println!();

    // let m_16: ExtendedMatrix<u16, f64> = ExtendedMatrix::create(3u16, 3u16,vec![3., -1., -2., 1., 7., -3., 3., -2., 10.]);
    // if let Ok((l_matrix, u_matrix)) = m_16.lu_decomposition()
    // {
    //     println!("{:?}", l_matrix.basic_matrix.extract_all_elements_values());
    //     l_matrix.show_matrix();
    //     println!();
    //     println!("{:?}", u_matrix.basic_matrix.extract_all_elements_values());
    //     u_matrix.show_matrix();
    //     println!();
    //     let m_19 = l_matrix.multiply_by_matrix(&u_matrix);
    //     match m_19
    //     {
    //         Ok(matrix) =>
    //             {
    //                 matrix.show_matrix();
    //                 // let basic_matrix = return_non_symmetric_matrix_struct(matrix.basic_matrix.clone());
    //                 // println!("{:?}\n", basic_matrix);
    //                 // println!("{:?}\n", mem::size_of_val(&basic_matrix));
    //             },
    //         Err(e) => println!("{}\n", e),
    //     }
    // }
    // println!();
    // if let Ok(determinant) = &m_16.determinant()
    // {
    //     println!("{}", determinant);
    //     println!();
    // }
    // if let Ok(inverse_matrix) = m_16.inverse()
    // {
    //     inverse_matrix.show_matrix();
    //     println!();
    //     if let Ok(unit_matrix) = m_16.multiply_by_matrix(&inverse_matrix)
    //     {
    //         unit_matrix.show_matrix();
    //     }
    // }
    //
    // println!();
    //
    // let x = 1.0;
    // let y = 1.0;
    // let z = 1.0;
    // let length = f64::sqrt(x * x + y * y + z * z);
    // let u = 0.0;
    // let v = 0.0;
    // let w = length;
    //
    // let alpha = f64::acos(((x * u + y * v + z * w) /
    //         (length * length)).into());
    //
    // let mut axis_x = y * length;
    // let mut axis_y = - x * length;
    // let mut axis_z = 0f64;
    // // if x != 0f64 && y == 0f64 && z == 0f64
    // // {
    // //     axis_z = x;
    // // }
    // // else
    // // {
    // //     axis_y = z * length;
    // //     axis_z = - y * length;
    // // }
    //
    // let norm = 1.0 / f64::sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
    //
    // let x_n = axis_x * norm;
    // let y_n = axis_y * norm;
    // let z_n = axis_z * norm;
    // let c = alpha.cos();
    // let s = alpha.sin();
    // let t = 1.0 - c;
    //
    // let q_11 = if (t * x_n * x_n + c).abs() < TOLERANCE { 0.0 } else { t * x_n * x_n + c };
    // let q_12 = if (t * x_n * y_n - z_n * s).abs() < TOLERANCE { 0.0 } else { t * x_n * y_n - z_n * s };
    // let q_13 = if (t * x_n * z_n + y_n * s).abs() < TOLERANCE { 0.0 } else { t * x_n * z_n + y_n * s };
    // let q_21 = if (t * x_n * y_n + z_n * s).abs() < TOLERANCE { 0.0 } else { t * x_n * y_n + z_n * s };
    // let q_22 = if (t * y_n * y_n + c).abs() < TOLERANCE { 0.0 } else { t * y_n * y_n + c };
    // let q_23 = if (t * y_n * z_n - x_n * s).abs() < TOLERANCE { 0.0 } else { t * y_n * z_n - x_n * s };
    // let q_31 = if (t * x_n * z_n - y_n * s).abs() < TOLERANCE { 0.0 } else { t * x_n * z_n - y_n * s };
    // let q_32 = if (t * y_n * z_n + x_n * s).abs() < TOLERANCE { 0.0 } else { t * y_n * z_n + x_n * s };
    // let q_33 = if (t * z_n * z_n + c).abs() < TOLERANCE { 0.0 } else { t * z_n * z_n + c };
    //
    // let v_1 = ExtendedMatrix::create(
    //     3u16, 1u16, vec![1., -1., 0.]);
    // let rotation_matrix = ExtendedMatrix::create(
    //     3u16, 3u16,
    //     vec![q_11, q_12, q_13, q_21, q_22, q_23, q_31, q_32, q_33]);
    //
    // println!("Rotation matrix:");
    // rotation_matrix.show_matrix();
    // println!();
    // if let Ok(m) = rotation_matrix.multiply_by_matrix(&v_1)
    // {
    //     println!("New p1:");
    //     m.show_matrix();
    //     println!();
    // }
    //
    // let v_2 = ExtendedMatrix::create(3u16, 1u16, vec![length, 0., 0.]);
    // if let Ok(inv_rot) = rotation_matrix.inverse()
    // {
    //     if let Ok(mut m) = inv_rot.multiply_by_matrix(&v_2)
    //     {
    //         m.show_matrix();
    //         m.transpose();
    //         m.show_matrix();
    //     }
    // }
    //
    // let mut m_20 = ExtendedMatrix::create(
    //     5u16, 5u16, vec![0., 0., 0., 0., 0.,
    //     0., 0., 0., 0., 0.,
    //     0., 1., 3., 0., 0.,
    //     0., 0., 2., 4., 0.,
    //     0., 0., 0., 0., 0.]);
    // // let mut m_20 = ExtendedMatrix::create(
    // //     5u16, 3u16, vec![0., 0., 0.,
    // //     0., 0., 0.,
    // //     0., 5., 3.,
    // //     0., 0., 2.,
    // //     0., 0., 0.]);
    // let zeros_rows_columns = m_20.remove_zeros_rows_columns();
    // println!("{:?}", zeros_rows_columns);
    // m_20.show_matrix();
    //
    //
    // println!();
    // // m_20.remove_zeros_row(2u16);
    //
    // m_20.show_matrix();
    // println!();
    //
    let node_1 = FeNode { number: 1u16, coordinates: GlobalCoordinates { x: 0.0, y: 0.0, z: 0.0 } };
    let node_2 = FeNode { number: 2u16, coordinates: GlobalCoordinates { x: 3.0, y: 3.0, z: 3.0 } };
    let node_3 = FeNode { number: 3u16, coordinates: GlobalCoordinates { x: 4.0, y: 5.0, z: 6.0 } };
    let mut elem_1 = Truss2n2ip::create(1u16, &node_2, &node_1, 1e6, 1.0, None);
    elem_1.state.rotation_matrix.show_matrix();
    println!();
    elem_1.update(&node_3, &node_2, 1e6, 1.0, None);
    elem_1.state.rotation_matrix.show_matrix();
    println!();
}
