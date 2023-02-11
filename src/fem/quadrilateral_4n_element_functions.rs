use extended_matrix::{FloatTrait, Vector3, VectorTrait, BasicOperationsTrait, Position};

use crate::fem::math_functions::compare_with_tolerance;


pub fn is_points_of_quadrilateral_on_the_same_line<V>(
    point_1: &[V], point_2: &[V], point_3: &[V], point_4: &[V],
) 
    -> bool
    where V: FloatTrait<Output = V>,
{
    let vector_1_2 = Vector3::create(
        &[point_2[0] - point_1[0], point_2[1] - point_1[1], point_2[2] - point_1[2]],
    );
    let vector_1_4 = Vector3::create(
        &[point_4[0] - point_1[0], point_4[1] - point_1[1], point_4[2] - point_1[2]],
    );

    let vector_2_1 = Vector3::create(
        &[point_1[0] - point_2[0], point_1[1] - point_2[1], point_1[2] - point_2[2]],
    );
    let vector_2_3 = Vector3::create(
        &[point_3[0] - point_2[0], point_3[1] - point_2[1], point_3[2] - point_2[2]],
    );

    let vector_3_2 = Vector3::create(
        &[point_2[0] - point_3[0], point_2[1] - point_3[1], point_2[2] - point_3[2]],
    );
    let vector_3_4 = Vector3::create(
        &[point_4[0] - point_3[0], point_4[1] - point_3[1], point_4[2] - point_3[2]],
    );

    let vector_4_3 = Vector3::create(
        &[point_3[0] - point_4[0], point_3[1] - point_4[1], point_3[2] - point_4[2]],
    );
    let vector_4_1 = Vector3::create(
        &[point_1[0] - point_4[0], point_1[1] - point_4[1], point_1[2] - point_4[2]],
    );

    let vector_pairs = [
        [vector_1_2, vector_1_4], [vector_2_1, vector_2_3], [vector_3_2, vector_3_4], [vector_4_3, vector_4_1],
    ];

    for vector_pair in vector_pairs
    {
        if vector_pair[0].cross_product(&vector_pair[1]).norm().expect("Norm could not be found!") == V::from(0f32)
        {
            return true;
        }
    }

    false
}


pub fn is_points_of_quadrilateral_on_the_same_plane<V>(
    point_1: &[V], point_2: &[V], point_3: &[V], point_4: &[V], abs_tol: V,
) 
    -> bool
    where V: FloatTrait<Output = V>,
{
    let vector_3_2 = Vector3::create(
        &[point_2[0] - point_3[0], point_2[1] - point_3[1], point_2[2] - point_3[2]],
    );
    let vector_3_4 = Vector3::create(
        &[point_4[0] - point_3[0], point_4[1] - point_3[1], point_4[2] - point_3[2]],
    ); 

    let normal_to_plane = vector_3_2.cross_product(&vector_3_4);
    let a = *normal_to_plane.get_element_value(&Position(0, 0)).expect("Incorrect position (0, 0)!");
    let b = *normal_to_plane.get_element_value(&Position(1, 0)).expect("Incorrect position (1, 0)!");
    let c = *normal_to_plane.get_element_value(&Position(2, 0)).expect("Incorrect position (0, 0)!");
    let d = V::from(-1f32) * (a * point_3[0] + b * point_3[1] + c * point_3[2]);
    
    if compare_with_tolerance(a * point_1[0] + b * point_1[1] + c * point_1[2] + d, abs_tol) != V::from(0f32)
    {
        return false
    }
    true
}
