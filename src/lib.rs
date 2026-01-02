//! finite_element_method
//!
//! A small Rust crate with finite-element method (FEM) building blocks used by the
//! `fea_app` demo project. It focuses on representing meshes/elements, assembling
//! stiffness matrices, and producing the data needed for solvers.
//!
//! This crate is intentionally **learning / exploration** oriented: the APIs are
//! kept straightforward and the code favors readability over maximum generality.
//!
//! ## Structure
//! - [`fem`] — high-level FEM entry point and helpers
//! - [`fem::structs`] — element implementations (truss/beam/plate), nodes, material props
//! - `convex_hull_on_plane` — small geometric helper used by some element routines
//!
//! The iterative/PCG solvers and the WebGPU compute path live in separate crates/repos.

mod fem;

pub use fem::{
    DOFParameter, ElementForceComponent, FEM, SeparatedStiffnessMatrix,
    SeparatedStiffnessMatrixSparse, convex_hull_on_four_points_on_plane,
    is_points_of_quadrilateral_on_the_same_line, is_points_of_quadrilateral_on_the_same_plane,
};

mod tests;
