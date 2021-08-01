pub trait FloatTrait
{
    fn powi(&self, n: i32) -> Self;
    fn sqrt(&self) -> Self;
    fn acos(&self) -> Self;
    fn cos(&self) -> Self;
    fn sin(&self) -> Self;
    fn abs(&self) -> Self;
}


impl FloatTrait for f32
{
    fn powi(&self, n: i32) -> Self
    {
        self.powi(n)
    }


    fn sqrt(&self) -> Self
    {
        self.sqrt()
    }


    fn acos(&self) -> Self
    {
        self.acos()
    }


    fn cos(&self) -> Self
    {
        self.cos()
    }


    fn sin(&self) -> Self
    {
        self.sin()
    }


    fn abs(&self) -> Self
    {
        self.abs()
    }
}


impl FloatTrait for f64
{
    fn powi(&self, n: i32) -> Self
    {
        self.powi(n)
    }


    fn sqrt(&self) -> Self
    {
        self.sqrt()
    }


    fn acos(&self) -> Self
    {
        self.acos()
    }


    fn cos(&self) -> Self
    {
        self.cos()
    }


    fn sin(&self) -> Self
    {
        self.sin()
    }


    fn abs(&self) -> Self
    {
        self.abs()
    }
}
