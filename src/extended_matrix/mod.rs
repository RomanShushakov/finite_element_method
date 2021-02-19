pub mod basic_matrix;
pub mod extended_matrix;
pub mod aux_functions_extended_matrix;

pub use crate::extended_matrix::extended_matrix::ExtendedMatrix;
pub use crate::extended_matrix::extended_matrix::Operation;
pub use crate::extended_matrix::aux_functions_extended_matrix::
    {
        matrices_dimensions_conformity_check, extract_element_value, remove_zero_values
    };
pub use crate::extended_matrix::basic_matrix::basic_matrix::BasicMatrix;
pub use crate::extended_matrix::basic_matrix::basic_matrix::
    {
        Shape, ZerosRowColumn, MatrixElementPosition
    };
pub use crate::extended_matrix::basic_matrix::basic_matrix::BasicMatrixType;
pub use crate::extended_matrix::basic_matrix::symmetric_matrix::SymmetricMatrix;
pub use crate::extended_matrix::basic_matrix::non_symmetric_matrix::NonSymmetricMatrix;
pub use crate::extended_matrix::basic_matrix::aux_functions_basic_matrix::
    {
        return_symmetric_matrix_struct, return_non_symmetric_matrix_struct, matrix_size_check,
        extract_value_by_index
    };

