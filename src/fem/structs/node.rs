use extended_matrix::FloatTrait;


pub const NODE_DOF: usize = 6;


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


    pub fn get_index(&self) -> usize
    {
        self.index
    }


    pub fn is_coordinates_same(&self, x: V, y: V, z: V) -> bool
    {
        self.x == x && self.y == y && self.z == z
    }


    pub fn get_coordinates(&self) -> [V; 3]
    {
        [self.x, self.y, self.z]
    }
}
