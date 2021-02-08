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


    let x = 0.0;
    let y = -1.0;
    let z = 1.0;
    let length = f64::sqrt(3. * 3. + 3. * 3. + 3. * 3.);
    let alpha = f64::acos(3. / length).to_degrees();

    let heading = f64::atan2(
        y * alpha.to_radians().sin() - x * z * (1.0 - alpha.to_radians().cos()),
        1.0 - (y * y + z * z ) * (1.0 - alpha.to_radians().cos()));
    println!("Heading: {}", heading.to_degrees());

    let attitude = f64::asin(
        x * y * (1.0 - alpha.to_radians().cos()) + z * alpha.to_radians().sin());
    println!("Attitude: {}", attitude.to_degrees());

    let bank = f64::atan2(
        x * alpha.to_radians().sin() - y * z * (1.0 - alpha.to_radians().cos()) ,
        1.0 - (x * x + z * z) * (1.0 - alpha.to_radians().cos()));
    println!("Bank: {}", bank.to_degrees());

    let q_11 = heading.cos() * attitude.cos();
    println!("Q_11: {}", q_11);
    let q_12 = - heading.cos() * attitude.sin() * bank.cos() + heading.sin() * bank.sin();
    println!("Q_12: {}", q_12);
    let q_13 = heading.cos() * attitude.sin() * bank.sin() + heading.sin() * bank.cos();
    println!("Q_13: {}", q_13);
    let q_21 = attitude.sin();
    println!("Q_21: {}", q_21);
    let q_22 = attitude.cos() * bank.cos();
    println!("Q_22: {}", q_22);
    let q_23 = - attitude.cos() * bank.sin();
    println!("Q_23: {}", q_23);
    let q_31 = - heading.sin() * attitude.cos();
    println!("Q_31: {}", q_31);
    let q_32 = heading.sin() * attitude.sin() * bank.cos() + heading.cos() * bank.sin();
    println!("Q_32: {}", q_32);
    let q_33 = - heading.sin() * attitude.sin() * bank.sin() + heading.cos() * bank.cos();
    println!("Q_33: {}", q_33);

    let v_1 = ExtendedMatrix::create(
        3u16, 1u16, vec![3., 3., 3.]);
    let v_2 = ExtendedMatrix::create(
        3u16, 1u16, vec![length, 0., 0.]);
    let rotational_matrix = ExtendedMatrix::create(
        3u16, 3u16,
        vec![q_11, q_12, q_13, q_21, q_22, q_23, q_31, q_32, q_33]);
    if let Ok(inv_m) = rotational_matrix.inverse()
    {
        if let Ok(m) = inv_m.multiply_by_matrix(&v_1)
        {
            m.show_matrix();
        }
    }

}
