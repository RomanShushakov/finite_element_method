use extended_matrix::SquareMatrix;


pub struct Truss<V>
{
    node_1_number: u32,
    node_2_number: u32,
    young_modulus: V,
    poisson_ratio: V,
    area: V,
    optional_area_2: Option<V>,
    rotation_matrix: SquareMatrix<V>,
    integration_points: [(V, V); 1],
    stiffness_matrix: SquareMatrix<V>,
}
