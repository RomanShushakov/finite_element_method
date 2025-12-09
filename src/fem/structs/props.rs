pub struct Props<V> {
    rel_tol: V,
    abs_tol: V,
    nodes_number: u32,
}

impl<V> Props<V>
where
    V: Copy,
{
    pub fn create(rel_tol: V, abs_tol: V, nodes_number: u32) -> Self {
        Props {
            rel_tol,
            abs_tol,
            nodes_number,
        }
    }

    pub fn get_rel_tol(&self) -> V {
        self.rel_tol
    }

    pub fn get_abs_tol(&self) -> V {
        self.abs_tol
    }

    pub fn get_nodes_number(&self) -> u32 {
        self.nodes_number
    }
}
