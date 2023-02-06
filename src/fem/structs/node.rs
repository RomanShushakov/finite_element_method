use extended_matrix::FloatTrait;


pub struct Node<V>
{
    index: usize,
    x: V,
    y: V,
    z: V,
}


impl<V> Node<V>
    where V: FloatTrait
{
    pub fn create(index: usize, x: V, y: V, z: V) -> Self
    {
        Node { index, x, y, z }
    }


    pub fn is_index_same(&self, index: usize) -> bool
    {
        self.index == index
    }


    pub fn is_coordinates_same(&self, x: V, y: V, z: V) -> bool
    {
        self.x == x && self.y == y && self.z == z
    }
}
