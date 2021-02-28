pub use crate::fem::element_analysis::fe_stress_strain_components::{StressStrainComponent};
pub use crate::fem::element_analysis::fe_force_moment_components::ForceComponent;
pub use crate::fem::element_analysis::fe_element_analysis_result::ElementAnalysisData;

pub use crate::fem::aux_functions_fe_model::compose_stiffness_sub_groups;
pub use crate::fem::fe_model::{FEModel, SeparatedMatrix};
pub use crate::fem::finite_elements::aux_functions_finite_elements::compare_with_tolerance;
pub use crate::fem::finite_elements::fe_node::{FeNode, GlobalCoordinates};
pub use crate::fem::finite_elements::finite_element::{FEData, FiniteElement};
pub use crate::fem::finite_elements::finite_element::FEType;
pub use crate::fem::finite_elements::finite_element::FiniteElementTrait;
pub use crate::fem::finite_elements::truss::truss2n2ip::Truss2n2ip;
pub use crate::fem::global_analysis::fe_boundary_condition::BCType;
pub use crate::fem::global_analysis::fe_boundary_condition::BoundaryCondition;
pub use crate::fem::global_analysis::fe_dof_parameter_data::{DOFParameterData, GlobalDOFParameter};
pub use crate::fem::global_analysis::fe_dof_parameter_data::GLOBAL_DOF;
pub use crate::fem::global_analysis::fe_global_analysis_result::GlobalAnalysisResult;
pub use crate::fem::global_analysis::fe_stiffness::STIFFNESS_TYPES_NUMBER;
pub use crate::fem::global_analysis::fe_stiffness::StiffnessGroup;
pub use crate::fem::global_analysis::fe_stiffness::StiffnessType;
pub use crate::fem::global_analysis::fe_global_analysis_result::Displacements;

pub mod finite_elements;
pub mod fe_model;
pub mod global_analysis;
pub mod element_analysis;

pub mod aux_functions_fe_model;

