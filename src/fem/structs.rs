//! Element and mesh data structures.
//!
//! This module contains the concrete element implementations (truss / beam / plate)
//! and small supporting types such as [`Node`] and [`Props`].
//!
//! The element types typically implement:
//! - local stiffness / transformation routines
//! - assembly helpers (mapping local DOFs to global DOFs)
//!
//! The goal is to keep each element self-contained so it can be tested/extended in
//! isolation.

mod props;
pub use props::Props;

mod node;
pub use node::{NODE_DOF, Node};

mod truss;
pub use truss::{TRUSS_NODE_DOF, Truss};

mod separated_stiffness_matrix;
pub use separated_stiffness_matrix::SeparatedStiffnessMatrix;

mod separated_stiffness_matrix_sparse;
pub use separated_stiffness_matrix_sparse::SeparatedStiffnessMatrixSparse;

mod beam;
pub use beam::{BEAM_NODE_DOF, Beam};

mod plate;
pub use plate::{PLATE_NODE_DOF, Plate};
