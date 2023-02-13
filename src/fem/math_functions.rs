use extended_matrix::FloatTrait;


pub fn compare_with_tolerance<V>(value: V, abs_tol: V) -> V
    where V: FloatTrait<Output = V>
{
    if value.my_abs() < abs_tol
    {
        V::from(0f32)
    }
    else
    {
        value
    }
}


pub fn power_func_x<V>(a: V, x: V, n: i32) -> V
    where V: FloatTrait<Output = V>
{
    (0..n).fold(a, |acc, _| acc * x)
}


pub fn derivative_x<V>(f: fn(V, V, i32) -> V, a: V, x: V, n: i32) -> V
    where V: FloatTrait<Output = V>
{
    let mut converted_n = V::from(0f32);
    (0..n).for_each(|_| converted_n += V::from(1f32));
    f(a * converted_n, x, n - 1)
}
