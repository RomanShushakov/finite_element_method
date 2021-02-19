pub mod finite_elements;
pub mod fe_model;

pub use crate::fem::finite_elements::fe::{StiffnessType, StiffnessGroup};
pub use crate::fem::finite_elements::fe_node::{FeNode, GlobalCoordinates};
pub use crate::fem::finite_elements::truss::truss2n2ip::Truss2n2ip;
pub use crate::fem::finite_elements::aux_functions_finite_elements::compare_with_tolerance;

pub use crate::fem::fe_model::FEModel;