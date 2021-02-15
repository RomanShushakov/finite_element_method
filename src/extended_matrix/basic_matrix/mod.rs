pub mod basic_matrix;
pub mod non_symmetric_matrix;
pub mod symmetric_matrix;
pub mod aux_functions_basic_matrix;


pub use basic_matrix::{BasicMatrix};
pub use basic_matrix::{Shape, MatrixElementPosition, ZerosRowColumn};
pub use basic_matrix::{BasicMatrixType};
pub use non_symmetric_matrix::NonSymmetricMatrix;
pub use symmetric_matrix::SymmetricMatrix;
pub use aux_functions_basic_matrix::
    {
        matrix_size_check, extract_value_by_index, return_symmetric_matrix_struct,
        return_non_symmetric_matrix_struct
    };
