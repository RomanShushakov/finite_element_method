#[derive(Debug)]
pub struct FENode<T, V>
{
    number: T,
    x: V,
    y: V,
    z: V,
}


impl<T, V> FENode<T, V>
    where T: PartialEq + Copy,
          V: PartialEq + Copy,
{
    pub fn create(number: T, x: V, y: V, z: V) -> Self
    {
        FENode { number, x, y, z }
    }


    pub fn update(&mut self, x: V, y: V, z: V)
    {
        self.x = x;
        self.y = y;
        self.z = z;
    }


    pub fn number_same(&self, number: T) -> bool
    {
        number == self.number
    }


    pub fn coordinates_same(&self, x: V, y: V, z: V) -> bool
    {
        (x, y, z) == (self.x, self.y, self.z)
    }


    pub fn extract_number(&self) -> T
    {
        self.number
    }


    pub fn extract_coordinates(&self) -> (V, V, V)
    {
        (self.x, self.y, self.z)
    }


    pub fn number(&self) -> T
    {
        self.number
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
