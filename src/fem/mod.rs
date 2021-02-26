pub mod finite_elements;
pub mod fe_model;
pub mod fe_boundary_condition;
pub mod aux_functions_fe_model;
pub mod fe_dof_parameter_data;
pub mod fe_stiffness;
pub mod fe_analysis_result;

pub use crate::fem::finite_elements::finite_element::{FiniteElementTrait};
pub use crate::fem::finite_elements::finite_element::{FEData, FiniteElement};
pub use crate::fem::finite_elements::finite_element::{FEType};
pub use crate::fem::finite_elements::fe_node::{FeNode, GlobalCoordinates};
pub use crate::fem::finite_elements::truss::truss2n2ip::Truss2n2ip;
pub use crate::fem::finite_elements::aux_functions_finite_elements::compare_with_tolerance;

pub use crate::fem::fe_model::{FEModel, SeparatedMatrix};

pub use crate::fem::fe_dof_parameter_data::{GlobalDOFParameter, DOFParameterData};
pub use crate::fem::fe_dof_parameter_data::{GLOBAL_DOF};

pub use crate::fem::fe_boundary_condition::{BoundaryCondition};
pub use crate::fem::fe_boundary_condition::{BCType};

pub use crate::fem::fe_stiffness::StiffnessGroup;
pub use crate::fem::fe_stiffness::StiffnessType;
pub use crate::fem::fe_stiffness::{STIFFNESS_TYPES_NUMBER};

pub use crate::fem::aux_functions_fe_model::compose_stiffness_sub_groups;

pub use crate::fem::fe_analysis_result::GlobalAnalysisResult;

