use crate::my_float::MyFloatTrait;


pub fn compare_with_tolerance<V>(value: V, tolerance: V) -> V
    where V: MyFloatTrait + PartialOrd + From<f32>
{
    if value.my_abs() < tolerance
    {
        V::from(0f32)
    }
    else
    {
        value
    }
}
