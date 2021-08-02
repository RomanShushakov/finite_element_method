use crate::float::MyFloatTrait;

use crate::{ElementsValues, TOLERANCE};


pub fn compare_with_tolerance<V>(value: V, tolerance: V) -> V
    where V: Default + MyFloatTrait + PartialOrd
{
    if value.my_abs() < tolerance
    {
        V::default()
    }
    else
    {
        value
    }
}
