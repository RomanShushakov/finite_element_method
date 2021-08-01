use crate::float::FloatTrait;

use crate::{ElementsValues, TOLERANCE};


pub fn compare_with_tolerance<V>(value: V, tolerance: V) -> V
    where V: Default + FloatTrait + PartialOrd
{
    if value.abs() < tolerance
    {
        V::default()
    }
    else
    {
        value
    }
}
