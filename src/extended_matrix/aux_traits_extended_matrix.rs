pub trait One
{
    fn one() -> Self;
}


impl One for f64
{
    fn one() -> f64
    {
        1.0
    }
}


impl One for f32
{
    fn one() -> f32
    {
        1.0
    }
}


impl One for i32
{
    fn one() -> i32
    {
        1
    }
}
