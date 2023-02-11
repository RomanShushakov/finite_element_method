pub mod finite_elements;
pub mod fe_model;
pub mod global_analysis;
pub mod element_analysis;
pub mod functions;
pub(super) mod separated_matrix;
pub(super) mod convex_hull_on_plane;

mod fem;
pub use fem::FEM;

mod structs;
pub use structs::SeparatedStiffnessMatrix;

mod methods_for_node_data_handle;

mod methods_for_truss_data_handle;

mod methods_for_bc_data_handle;
pub use methods_for_bc_data_handle::DOFParameter;

mod methods_for_separate_stiffness_matrix;

mod methods_for_global_analysis;

mod methods_for_element_analysis;

mod methods_for_beam_data_handle;

mod math_functions;

mod bar_2n_element_functions;

mod quadrilateral_4n_element_functions;
pub use quadrilateral_4n_element_functions::
{
    is_points_of_quadrilateral_on_the_same_line, is_points_of_quadrilateral_on_the_same_plane,
};
