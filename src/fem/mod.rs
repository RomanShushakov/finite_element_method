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
pub use methods_for_element_analysis::ElementForceComponent;

mod methods_for_beam_data_handle;

mod math_functions;

mod bar_2n_element_functions;

mod convex_hull_on_plane;

mod quadrilateral_4n_element_functions;
pub use quadrilateral_4n_element_functions::
{
    is_points_of_quadrilateral_on_the_same_line, is_points_of_quadrilateral_on_the_same_plane,
    find_rotation_matrix_elements_of_quadrilateral, convex_hull_on_four_points_on_plane,
};

mod methods_for_plate_data_handle;
