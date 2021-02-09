pub struct GlobalCoordinates<T>
{
    pub x: T,
    pub y: T,
    pub z: T,
}


pub struct FeNode<T, V>
{
    pub number: T,
    pub coordinates: GlobalCoordinates<V>
}
