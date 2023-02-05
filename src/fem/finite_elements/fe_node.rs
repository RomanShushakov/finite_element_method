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
    pub(crate) fn create(x: V, y: V, z: V) -> Self
    {
        FENode { x, y, z }
    }


    pub(crate) fn update(&mut self, x: V, y: V, z: V)
    {
        self.x = x;
        self.y = y;
        self.z = z;
    }


    pub(crate) fn is_coordinates_same(&self, x: V, y: V, z: V) -> bool
    {
        (x, y, z) == (self.x, self.y, self.z)
    }

    
    pub(crate) fn copy_coordinates(&self) -> (V, V, V)
    {
        (self.x, self.y, self.z)
    }


    pub(crate) fn copy_x(&self) -> V
    {
        self.x
    }


    pub(crate) fn copy_y(&self) -> V
    {
        self.y
    }


    pub(crate) fn copy_z(&self) -> V
    {
        self.z
    }
}


#[derive(Debug, Clone)]
pub struct DeletedFENodeData<V>
{
    number: u32,
    x: V,
    y: V,
    z: V,
}


impl<V> DeletedFENodeData<V>
    where V: Copy + PartialEq,
{
    pub(crate) fn create(number: u32, deleted_node: FENode<V>) -> Self
    {
        let (x, y, z) = deleted_node.copy_coordinates();
        DeletedFENodeData { number, x, y, z }
    }


    pub fn copy_number(&self) -> u32
    {
        self.number
    }


    pub fn copy_coordinates(&self) -> (V, V, V)
    {
        (self.x, self.y, self.z)
    }
}
