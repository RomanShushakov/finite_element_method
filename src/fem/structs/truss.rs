use std::fmt::Debug;

use extended_matrix::{SquareMatrix, FloatTrait};


enum TrussDataError<V>
{
    YoungModulus(V),
    PoissonRatio(V),
    Area(V),
    Area2(V),
}


impl<V> TrussDataError<V>
    where V: Debug
{
    fn compose_error_message(&self) -> String
    {
        match self
        {
            TrussDataError::YoungModulus(value) => format!("Young's modulus {value:?} is less or equal to zero!"),
            TrussDataError::PoissonRatio(value) => format!("Poisson's ratio {value:?} is less or equal to zero!"),
            TrussDataError::Area(value) => format!("Area {value:?} is less or equal to zero!"),
            TrussDataError::Area2(value) => format!("Area2 {value:?} is less or equal to zero!"),
        }
    }
}


fn check_truss_properties<V>(
    young_modulus: V, poisson_ratio: V, area: V, optional_area_2: Option<V>,
) 
    -> Result<(), String>
    where V: Debug + PartialEq + PartialOrd + From<f32>
{
    if young_modulus <= V::from(0f32)
    {
        return Err(TrussDataError::<V>::YoungModulus(young_modulus).compose_error_message());
    }
    if poisson_ratio <= V::from(0f32)
    {
        return Err(TrussDataError::<V>::PoissonRatio(poisson_ratio).compose_error_message());
    }
    if area <= V::from(0f32)
    {
        return Err(TrussDataError::<V>::Area(area).compose_error_message());
    }
    if let Some(area_2) = optional_area_2
    {
        if area_2 <= V::from(0f32)
        {
            return Err(TrussDataError::<V>::Area2(area_2).compose_error_message());
        }
    }
    Ok(())
}


pub struct Truss<V>
{
    node_1_number: u32,
    node_2_number: u32,
    young_modulus: V,
    poisson_ratio: V,
    area: V,
    optional_area_2: Option<V>,
    rotation_matrix_elements: [V; 9],
    integration_points: [(V, V); 1],
}


impl<V> Truss<V>
{
    pub fn create(
        node_1_number: u32,
        node_2_number: u32,
        young_modulus: V,
        poisson_ratio: V,
        area: V,
        optional_area_2: Option<V>,
    )
        -> Result<Self, String>
        where V: FloatTrait
    {
        check_truss_properties(young_modulus, poisson_ratio, area, optional_area_2)?;

        let rotation_matrix_elements = [V::from(0f32); 9];
        let integration_points = [(V::from(0f32), V::from(2f32))];

        Ok(
            Truss
            {
                node_1_number, node_2_number, young_modulus, poisson_ratio, area, optional_area_2, 
                rotation_matrix_elements, integration_points,
            }
        )
    }
}
