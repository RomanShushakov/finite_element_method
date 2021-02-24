pub mod finite_elements;
pub mod fe_model;
pub mod fe_force_displacement;
pub mod aux_functions_fe_model;

pub use crate::fem::finite_elements::finite_element::{FiniteElementTrait};
pub use crate::fem::finite_elements::finite_element::{StiffnessGroup, FEData, FiniteElement};
pub use crate::fem::finite_elements::finite_element::{StiffnessType, FEType};
pub use crate::fem::finite_elements::finite_element::{STIFFNESS_TYPES_NUMBER};
pub use crate::fem::finite_elements::fe_node::{FeNode, GlobalCoordinates};
pub use crate::fem::finite_elements::truss::truss2n2ip::Truss2n2ip;
pub use crate::fem::finite_elements::aux_functions_finite_elements::compare_with_tolerance;

pub use crate::fem::fe_model::FEModel;
pub use crate::fem::fe_model::GLOBAL_DOF;

pub use crate::fem::fe_force_displacement::{DOFParameterData, Force, Displacement};
pub use crate::fem::fe_force_displacement::{GlobalDOFParameter};

pub use crate::fem::aux_functions_fe_model::compose_stiffness_sub_groups;

