#[derive(Debug)]
pub struct GlobalCoordinates<T>
{
    pub x: T,
    pub y: T,
    pub z: T,
}


#[derive(Debug)]
pub struct FeNode<T, V>
{
    pub number: T,
    pub coordinates: GlobalCoordinates<V>
}


impl<T, V> FeNode<T, V>
{
    pub fn create(number: T, x: V, y: V, z: V) -> Self
    {
        FeNode { number, coordinates: GlobalCoordinates { x, y, z } }
    }


    pub fn update(&mut self, x: V, y: V, z: V)
    {
        self.coordinates = GlobalCoordinates { x, y, z };
    }
}
