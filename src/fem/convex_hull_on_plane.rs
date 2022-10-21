use std::fmt::Debug;
use std::cmp::Ordering;

use extended_matrix::traits::{UIntTrait, FloatTrait};


pub fn quick_sort<T>(arr: &mut [T])
    where T: PartialOrd
{
    let len = arr.len();
    _quick_sort(arr, 0, (len - 1) as isize);
}


fn _quick_sort<T>(arr: &mut [T], low: isize, high: isize)
    where T: PartialOrd
{
    if low < high
    {
        let p = partition(arr, low, high);
        _quick_sort(arr, low, p - 1);
        _quick_sort(arr, p + 1, high);
    }
}


fn partition<T>(arr: &mut [T], low: isize, high: isize) -> isize
    where T: PartialOrd
{
    let pivot = high as usize;
    let mut store_index = low - 1;
    let mut last_index = high;

    loop
    {
        store_index += 1;
        while arr[store_index as usize] < arr[pivot]
        {
            store_index += 1;
        }
        last_index -= 1;
        while last_index >= 0 && arr[last_index as usize] > arr[pivot]
        {
            last_index -= 1;
        }
        if store_index >= last_index
        {
            break;
        }
        else
        {
            arr.swap(store_index as usize, last_index as usize);
        }
    }
    arr.swap(store_index as usize, pivot as usize);
    store_index
}


#[derive(Debug, Clone)]
pub struct Point<T, V>
{
    n: V,
    x: T,
    y: T,
}


impl<T, V> Point<T, V>
    where V: UIntTrait<Output = V>,
          T: FloatTrait<Output = T, Other = T>
{
    pub fn create(n: V, x: T, y: T) -> Self
    {
        Point { n, x, y }
    }


    fn vector_from_origin_to_point(&self) -> Vector<T, V>
    {
        Vector::create(
            V::from(0u8), T::from(0.0), T::from(0.0),
            self.n, self.x, self.y
        )
    }


    pub fn copy_number(&self) -> V
    {
        self.n
    }
}


impl<T, V> PartialOrd for Point<T, V>
    where V: UIntTrait<Output = V>,
          T: FloatTrait<Output = T, Other = T>
{
    fn partial_cmp(&self, other: &Point<T, V>) -> Option<Ordering>
    {
        let directional_vector: Vector<T, V> = Vector::create_directional_vector();
        let lhs_vector_from_origin = self.vector_from_origin_to_point();
        let rhs_vector_from_origin = other.vector_from_origin_to_point();
        directional_vector.cosine_of_angle_between_vectors(&lhs_vector_from_origin).my_acos().my_to_degrees()
            .partial_cmp(&(directional_vector.cosine_of_angle_between_vectors(&rhs_vector_from_origin).my_acos())
                .my_to_degrees())
    }
}


impl<T, V> PartialEq for Point<T, V>
    where V: UIntTrait<Output = V>,
          T: FloatTrait<Output = T, Other = T>
{
    fn eq(&self, other: &Point<T, V>) -> bool
    {
        let directional_vector = Vector::create_directional_vector();
        let lhs_vector_from_origin = self.vector_from_origin_to_point();
        let rhs_vector_from_origin = other.vector_from_origin_to_point();
        directional_vector.cosine_of_angle_between_vectors(&lhs_vector_from_origin).my_acos().my_to_degrees() ==
            directional_vector.cosine_of_angle_between_vectors(&rhs_vector_from_origin).my_acos().my_to_degrees()

    }
}


#[derive(Debug)]
struct Vector<T, V>
{
    point_1: Point<T, V>,
    point_2: Point<T, V>,
}


fn double_signed_area<T, V>(p_1: &Point<T, V>, p_2: &Point<T, V>, p_3: &Point<T, V>) -> T
    where T: FloatTrait<Output = T, Other = T>
{
    (p_2.x - p_1.x) * (p_3.y - p_1.y) - (p_2.y - p_1.y) * (p_3.x - p_1.x)
}


impl<T, V> Vector<T, V>
    where V: UIntTrait<Output = V>,
          T: FloatTrait<Output = T, Other = T>
{
    fn create_directional_vector() -> Self
    {
        Vector::create(
            V::from(0u8), T::from(0.0), T::from(0.0),
            V::from(0u8), T::from(1.0), T::from(0.0)
        )
    }


    fn create(n_1: V, x_1: T, y_1: T, n_2: V, x_2: T, y_2: T) -> Self
    {
        let point_1 = Point { n: n_1, x: x_1, y: y_1 };
        let point_2 = Point { n: n_2, x: x_2, y: y_2 };
        Vector { point_1, point_2 }
    }

    fn vector_length(&self) -> T
    {
        ((self.point_1.x - self.point_2.x) * (self.point_1.x - self.point_2.x) +
            (self.point_1.y - self.point_2.y) * (self.point_1.y - self.point_2.y)).my_sqrt()
    }


    fn scalar_product(&self, other: &Vector<T, V>) -> T
    {
        (self.point_1.x - self.point_2.x) * (other.point_1.x - other.point_2.x) +
        (self.point_1.y - self.point_2.y) * (other.point_1.y - other.point_2.y)
    }


    fn cosine_of_angle_between_vectors(&self, other: &Vector<T, V>) -> T
    {
        let lhs_length = self.vector_length();
        let rhs_length = other.vector_length();
        let scalar_product = self.scalar_product(other);
        scalar_product / (lhs_length * rhs_length)
    }
}


pub fn convex_hull_on_plane<T, V>(data: &[Point<T, V>]) -> Vec<Point<T, V>>
    where V: UIntTrait<Output = V>,
          T: FloatTrait<Output = T, Other = T>
{
    let mut updated_data = data.to_vec();
    let mut shift_x = updated_data[0].x;
    let mut min_y = updated_data[0].y;
    let mut min_y_position = 0;
    for i in 0..updated_data.len()
    {
        if updated_data[i].y < min_y
        {
            shift_x = updated_data[i].x;
            min_y = updated_data[i].y;
            min_y_position = i;
        }
    }
    updated_data.swap(0, min_y_position);
    for i in 0..updated_data.len()
    {
        updated_data[i].x -= shift_x;
        updated_data[i].y -= min_y;
    }
    quick_sort(&mut updated_data[1..]);
    let mut i = 0;
    while i + 2 < updated_data.len()
    {
        let doubled_area = double_signed_area(
            &updated_data[i], &updated_data[i + 1], &updated_data[i + 2]);
        if doubled_area <= T::from(0.0)
        {
            updated_data.remove(i + 1);
        }
        else
        {
            i += 1;
        }
    }
    for i in 0..updated_data.len()
    {
        updated_data[i].x += shift_x;
        updated_data[i].y += min_y;
    }
    updated_data
}
