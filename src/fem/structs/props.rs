pub struct Props<V>
{
    rel_tol: V,
    abs_tol: V,
    nodes_number: u32,
}


impl<V> Props<V>
{
    pub fn create(rel_tol: V, abs_tol: V, nodes_number: u32) -> Self
    {
        Props { rel_tol, abs_tol, nodes_number }
    }
}
