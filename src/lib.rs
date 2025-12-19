mod fem;

pub use fem::{
    DOFParameter, ElementForceComponent, FEM, SeparatedStiffnessMatrix,
    SeparatedStiffnessMatrixSparse, convex_hull_on_four_points_on_plane,
    is_points_of_quadrilateral_on_the_same_line, is_points_of_quadrilateral_on_the_same_plane,
};

mod tests;
