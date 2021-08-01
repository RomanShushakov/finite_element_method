use crate::{ElementsNumbers, ElementsValues};

#[derive(Debug, PartialEq)]
pub struct GlobalCoordinates<V>
{
    pub x: V,
    pub y: V,
    pub z: V,
}


impl<V> GlobalCoordinates<V>
    where V: Copy
{
    fn extract(&self) -> (V, V, V)
    {
        (self.x, self.y, self.z)
    }
}


#[derive(Debug)]
pub struct FENode<T, V>
{
    pub number: T,
    pub coordinates: GlobalCoordinates<V>
}


impl<T, V> FENode<T, V>
    where T: PartialEq + Copy,
          V: PartialEq + Copy,
{
    pub fn create(number: T, x: V, y: V, z: V) -> Self
    {
        FENode { number, coordinates: GlobalCoordinates { x, y, z } }
    }


    pub fn update(&mut self, x: V, y: V, z: V)
    {
        self.coordinates = GlobalCoordinates { x, y, z };
    }


    pub fn number_same(&self, number: T) -> bool
    {
        number == self.number
    }


    pub fn coordinates_same(&self, x: V, y: V, z: V) -> bool
    {
        GlobalCoordinates { x, y, z } == self.coordinates
    }


    pub fn extract_number(&self) -> T
    {
        self.number
    }


    pub fn extract_coordinates(&self) -> (V, V, V)
    {
        self.coordinates.extract()
    }

}
