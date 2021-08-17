pub trait MyFloatTrait
{
    type Other;

    fn my_powi(&self, n: i32) -> Self;
    fn my_sqrt(&self) -> Self;
    fn my_acos(&self) -> Self;
    fn my_cos(&self) -> Self;
    fn my_sin(&self) -> Self;
    fn my_abs(&self) -> Self;
    fn my_asin(&self) -> Self;
    fn my_atan2(&self, other: Self::Other) -> Self;
    fn my_atan(&self) -> Self;
    fn my_to_degrees(&self) -> Self;
}


impl MyFloatTrait for f32
{
    type Other = f32;

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


    fn my_asin(&self) -> Self
    {
        self.asin()
    }


    fn my_atan2(&self, other: Self::Other) -> Self
    {
        self.atan2(other)
    }


    fn my_atan(&self) -> Self
    {
        self.atan()
    }


    fn my_to_degrees(&self) -> Self
    {
        self.to_degrees()
    }
}


impl MyFloatTrait for f64
{
    type Other = f64;

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


    fn my_asin(&self) -> Self
    {
        self.asin()
    }


    fn my_atan2(&self, other: Self::Other) -> Self
    {
        self.atan2(other)
    }


    fn my_atan(&self) -> Self
    {
        self.atan()
    }


    fn my_to_degrees(&self) -> Self
    {
        self.to_degrees()
    }
}
