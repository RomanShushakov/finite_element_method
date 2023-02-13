mod fem;

pub use fem::
{
    FEM, DOFParameter, SeparatedStiffnessMatrix, is_points_of_quadrilateral_on_the_same_line,
    is_points_of_quadrilateral_on_the_same_plane, convex_hull_on_four_points_on_plane,
    ElementForceComponent,
};
