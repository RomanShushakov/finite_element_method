use crate::{ElementsNumbers, ElementsValues};

#[derive(Debug, PartialEq)]
pub struct GlobalCoordinates<T>
{
    pub x: T,
    pub y: T,
    pub z: T,
}


impl<T> GlobalCoordinates<T>
    where T: Copy + Into<ElementsValues>
{
    fn extract(&self) -> (ElementsValues, ElementsValues, ElementsValues)
    {
        (self.x.into(), self.y.into(), self.z.into())
    }
}


#[derive(Debug)]
pub struct FENode<T, V>
{
    pub number: T,
    pub coordinates: GlobalCoordinates<V>
}


impl<T, V> FENode<T, V>
    where T: PartialEq + Copy + Into<ElementsNumbers>,
          V: PartialEq + Copy + Into<ElementsValues>,
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


    pub fn extract_number(&self) -> ElementsNumbers
    {
        self.number.into()
    }


    pub fn extract_coordinates(&self) -> (ElementsValues, ElementsValues, ElementsValues)
    {
        self.coordinates.extract()
    }

}
