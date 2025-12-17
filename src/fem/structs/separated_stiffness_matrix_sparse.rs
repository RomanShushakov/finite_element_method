#[derive(Clone, Debug)]
pub struct SeparatedStiffnessMatrixSparse<V> {
    k_aa_indexes: Vec<usize>,
    k_bb_indexes: Vec<usize>,
    k_aa_triplets: Vec<(usize, usize, V)>,
    k_ab_triplets: Vec<(usize, usize, V)>,
    k_ba_triplets: Vec<(usize, usize, V)>,
    k_bb_triplets: Vec<(usize, usize, V)>,
    n_aa: usize,
    n_bb: usize,
}

impl<V> SeparatedStiffnessMatrixSparse<V> {
    pub fn create(
        k_aa_indexes: Vec<usize>,
        k_bb_indexes: Vec<usize>,
        k_aa_triplets: Vec<(usize, usize, V)>,
        k_ab_triplets: Vec<(usize, usize, V)>,
        k_ba_triplets: Vec<(usize, usize, V)>,
        k_bb_triplets: Vec<(usize, usize, V)>,
        n_aa: usize,
        n_bb: usize,
    ) -> Self {
        SeparatedStiffnessMatrixSparse {
            k_aa_indexes,
            k_bb_indexes,
            k_aa_triplets,
            k_ab_triplets,
            k_ba_triplets,
            k_bb_triplets,
            n_aa,
            n_bb,
        }
    }

    pub fn get_k_aa_indexes(&self) -> &Vec<usize> {
        &self.k_aa_indexes
    }

    pub fn get_k_bb_indexes(&self) -> &Vec<usize> {
        &self.k_bb_indexes
    }

    pub fn get_k_aa_triplets(&self) -> &Vec<(usize, usize, V)> {
        &self.k_aa_triplets
    }

    pub fn get_k_ab_triplets(&self) -> &Vec<(usize, usize, V)> {
        &self.k_ab_triplets
    }

    pub fn get_k_ba_triplets(&self) -> &Vec<(usize, usize, V)> {
        &self.k_ba_triplets
    }

    pub fn get_k_bb_triplets(&self) -> &Vec<(usize, usize, V)> {
        &self.k_bb_triplets
    }

    pub fn get_n_aa(&self) -> usize {
        self.n_aa
    }

    pub fn get_n_bb(&self) -> usize {
        self.n_bb
    }
}
