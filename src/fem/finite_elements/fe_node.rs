#[derive(Debug)]
pub struct FENode<V>
{
    x: V,
    y: V,
    z: V,
}


impl<V> FENode<V>
    where V: PartialEq + Copy,
{
    pub fn create(x: V, y: V, z: V) -> Self
    {
        FENode { x, y, z }
    }


    pub fn update(&mut self, x: V, y: V, z: V)
    {
        self.x = x;
        self.y = y;
        self.z = z;
    }


    pub fn is_coordinates_same(&self, x: V, y: V, z: V) -> bool
    {
        (x, y, z) == (self.x, self.y, self.z)
    }

    
    pub fn extract_coordinates(&self) -> (V, V, V)
    {
        (self.x, self.y, self.z)
    }


    pub fn x(&self) -> V
    {
        self.x
    }


    pub fn y(&self) -> V
    {
        self.y
    }


    pub fn z(&self) -> V
    {
        self.z
    }
}
