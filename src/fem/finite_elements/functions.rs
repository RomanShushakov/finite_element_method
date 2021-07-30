use crate::{ElementsValues, TOLERANCE};


pub fn compare_with_tolerance(value: ElementsValues) -> ElementsValues
{
    if value.abs() < TOLERANCE { 0.0 } else { value }
}