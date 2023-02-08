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
pub use structs::Props;

mod methods_for_node_data_handle;

mod methods_for_truss_data_handle;

mod methods_for_bc_data_handle;
pub use methods_for_bc_data_handle::DOFParameter;

mod methods_for_shrink_indexes;
