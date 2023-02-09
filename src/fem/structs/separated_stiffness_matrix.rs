use extended_matrix::{Matrix, SquareMatrix};


pub struct SeparatedStiffnessMatrix<V>
{
    k_aa_indexes: Vec<usize>,
    k_bb_indexes: Vec<usize>,
    k_aa_skyline: Vec<usize>,
    k_aa_matrix: SquareMatrix<V>,
    k_ab_matrix: Matrix<V>,
    k_ba_matrix: Matrix<V>,
    k_bb_matrix: SquareMatrix<V>,
}


impl<V> SeparatedStiffnessMatrix<V>
{
    pub fn create(
        k_aa_indexes: Vec<usize>,
        k_bb_indexes: Vec<usize>,
        k_aa_skyline: Vec<usize>,
        k_aa_matrix: SquareMatrix<V>,
        k_ab_matrix: Matrix<V>,
        k_ba_matrix: Matrix<V>,
        k_bb_matrix: SquareMatrix<V>,
    )
        -> Self
    {
        SeparatedStiffnessMatrix 
        { 
            k_aa_indexes,
            k_bb_indexes,
            k_aa_skyline,
            k_aa_matrix,
            k_ab_matrix,
            k_ba_matrix,
            k_bb_matrix
        }
    }


    pub fn get_k_aa_indexes(&self) -> &Vec<usize>
    {
        &self.k_aa_indexes
    }


    pub fn get_k_bb_indexes(&self) -> &Vec<usize>
    {
        &self.k_bb_indexes
    }


    pub fn get_k_aa_skyline(&self) -> &Vec<usize>
    {
        &self.k_aa_skyline
    }


    pub fn get_k_aa_matrix(&self) -> &SquareMatrix<V>
    {
        &self.k_aa_matrix
    }


    pub fn get_k_ab_matrix(&self) -> &Matrix<V>
    {
        &self.k_ab_matrix
    }


    pub fn get_k_ba_matrix(&self) -> &Matrix<V>
    {
        &self.k_ba_matrix
    }


    pub fn get_k_bb_matrix(&self) -> &SquareMatrix<V>
    {
        &self.k_bb_matrix
    }
}
