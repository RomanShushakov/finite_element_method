use extended_matrix::extended_matrix::ExtendedMatrix;

pub struct SeparatedMatrix<T, V>
{
    k_aa: ExtendedMatrix<T, V>,
    k_ab: ExtendedMatrix<T, V>,
    k_ba: ExtendedMatrix<T, V>,
    k_bb: ExtendedMatrix<T, V>,
}


impl<T, V> SeparatedMatrix<T, V>
{
    pub fn create(k_aa: ExtendedMatrix<T, V>, k_ab: ExtendedMatrix<T, V>,
        k_ba: ExtendedMatrix<T, V>, k_bb: ExtendedMatrix<T, V>) -> Self
    {
        SeparatedMatrix { k_aa, k_ab, k_ba, k_bb }
    }


    pub fn k_aa(&self) -> &ExtendedMatrix<T, V>
    {
        &self.k_aa
    }


    pub fn k_ab(&self) -> &ExtendedMatrix<T, V>
    {
        &self.k_ab
    }


    pub fn k_ba(&self) -> &ExtendedMatrix<T, V>
    {
        &self.k_ba
    }


    pub fn k_bb(&self) -> &ExtendedMatrix<T, V>
    {
        &self.k_bb
    }
}
