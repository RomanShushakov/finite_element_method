pub trait MyFloatTrait
{
    fn my_powi(&self, n: i32) -> Self;
    fn my_sqrt(&self) -> Self;
    fn my_acos(&self) -> Self;
    fn my_cos(&self) -> Self;
    fn my_sin(&self) -> Self;
    fn my_abs(&self) -> Self;
}


impl MyFloatTrait for f32
{
    fn my_powi(&self, n: i32) -> Self
    {
        self.powi(n)
    }


    fn my_sqrt(&self) -> Self
    {
        self.sqrt()
    }


    fn my_acos(&self) -> Self
    {
        self.acos()
    }


    fn my_cos(&self) -> Self
    {
        self.cos()
    }


    fn my_sin(&self) -> Self
    {
        self.sin()
    }


    fn my_abs(&self) -> Self
    {
        self.abs()
    }
}


impl MyFloatTrait for f64
{
    fn my_powi(&self, n: i32) -> Self
    {
        self.powi(n)
    }


    fn my_sqrt(&self) -> Self
    {
        self.sqrt()
    }


    fn my_acos(&self) -> Self
    {
        self.acos()
    }


    fn my_cos(&self) -> Self
    {
        self.cos()
    }


    fn my_sin(&self) -> Self
    {
        self.sin()
    }


    fn my_abs(&self) -> Self
    {
        self.abs()
    }
}
