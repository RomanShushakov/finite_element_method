use extended_matrix::{Matrix, SquareMatrix};

pub struct SeparatedMatrix<V>
{
    k_aa: SquareMatrix<V>,
    k_ab: Matrix<V>,
    k_ba: Matrix<V>,
    k_bb: SquareMatrix<V>,
}


impl<V> SeparatedMatrix<V>
{
    pub fn create(
        k_aa: SquareMatrix<V>, 
        k_ab: Matrix<V>,
        k_ba: Matrix<V>, 
        k_bb: SquareMatrix<V>
    ) 
        -> Self
    {
        SeparatedMatrix { k_aa, k_ab, k_ba, k_bb }
    }


    pub fn ref_k_aa(&self) -> &SquareMatrix<V>
    {
        &self.k_aa
    }


    pub fn ref_k_ab(&self) -> &Matrix<V>
    {
        &self.k_ab
    }


    pub fn ref_k_ba(&self) -> &Matrix<V>
    {
        &self.k_ba
    }


    pub fn ref_k_bb(&self) -> &SquareMatrix<V>
    {
        &self.k_bb
    }
}
