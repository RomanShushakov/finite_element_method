pub trait MinusOne
{
    fn minus_one() -> Self;
}


impl MinusOne for f32
{
    fn minus_one() -> Self
    {
        -1f32
    }
}


impl MinusOne for f64
{
    fn minus_one() -> Self
    {
        -1f64
    }
}
