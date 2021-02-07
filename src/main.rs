mod extended_matrix;
use extended_matrix::extended_matrix::{ExtendedMatrix};
use extended_matrix::basic_matrix::return_symmetric_matrix_struct;

use std::mem;
use crate::extended_matrix::basic_matrix::return_non_symmetric_matrix_struct;


pub type ElementsNumbers = u16;


fn main()
{
    // let m_10 = ExtendedMatrix::create(3u16, 3u16, vec![1.0, 0.0, 0.0, 0.03333333333333333, 1.0, 0.0, 0.09999999999999999, -0.027129938124702525, 1.0]);
    // println!("{:?}\n", m_10.basic_matrix.define_type());
    // let m_11 = ExtendedMatrix::create(3u16, 3u16, vec![3.0, -0.1, -0.2, 0.0, 7.003333333333333, -0.29333333333333333, 0.0, 0.0, 10.012041884816753]);
    // let m_12 = &m_10.multiply_by_matrix(&m_11);
    // match m_12
    // {
    //     Ok(matrix) =>
    //         {
    //             matrix.show_matrix();
    //             // let basic_matrix = return_symmetric_matrix_struct(matrix.basic_matrix.clone());
    //             // println!("{:?}\n", basic_matrix);
    //             // println!("{:?}\n", mem::size_of_val(&basic_matrix));
    //         },
    //     Err(e) => println!("{}\n", e),
    // }
    // println!();
    //
    // let m_13 = ExtendedMatrix::create(
    //     3u16, 3u16,
    //     vec![3.0, -0.1, -0.2, 0.1, 7.0, -0.3, 0.3, -0.2, 10.0]);
    //
    // let m_14 = ExtendedMatrix::create(
    //     3u16, 1u16, vec![7.85, -19.3, 71.4]);
    // if let Ok(m_15) = &m_13.naive_gauss_elimination(&m_14)
    // {
    //     m_15.show_matrix();
    // }
    // println!();

    let m_16: ExtendedMatrix<u16, f64> = ExtendedMatrix::create(3u16, 3u16,vec![3., -1., -2., 1., 7., -3., 3., -2., 10.]);
    if let Ok((l_matrix, u_matrix)) = m_16.lu_decomposition()
    {
        println!("{:?}", l_matrix.basic_matrix.extract_all_elements_values());
        l_matrix.show_matrix();
        println!();
        println!("{:?}", u_matrix.basic_matrix.extract_all_elements_values());
        u_matrix.show_matrix();
        println!();
        let m_19 = l_matrix.multiply_by_matrix(&u_matrix);
        match m_19
        {
            Ok(matrix) =>
                {
                    matrix.show_matrix();
                    // let basic_matrix = return_non_symmetric_matrix_struct(matrix.basic_matrix.clone());
                    // println!("{:?}\n", basic_matrix);
                    // println!("{:?}\n", mem::size_of_val(&basic_matrix));
                },
            Err(e) => println!("{}\n", e),
        }
    }
    println!();
    if let Ok(determinant) = &m_16.determinant()
    {
        println!("{}", determinant);
        println!();
    }
    if let Ok(inverse_matrix) = m_16.inverse()
    {
        inverse_matrix.show_matrix();
        println!();
        if let Ok(unit_matrix) = m_16.multiply_by_matrix(&inverse_matrix)
        {
            unit_matrix.show_matrix();
        }
    }
}
