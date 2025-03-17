use crate::ecc::{EllipticCurve};
use crate::point::Point;
use num_traits::{FromPrimitive, One, ToPrimitive, Zero};
use num_integer::Integer;
use rand::Rng;
use std::ops::{Add, Sub, Mul, Rem, Div};

pub struct DSA<T> {
    pub g: Point<T>,
    pub n: T,
    pub ec: Box<dyn EllipticCurve<T>>,
}

impl<T> DSA<T>
where
    T: One
    + Zero
    + Copy
    + PartialEq
    + Clone
    + PartialOrd
    + FromPrimitive
    + ToPrimitive
    + Integer
    + From<u8>
    + Add<Output = T>
    + Sub<Output = T>
    + Mul<Output = T>
    + Rem<Output = T>

    + Div<Output = T>,
{
    /// Creates a new DSA instance.
    ///
    /// It verifies that the generator `g` is a valid point on the curve `ec`
    /// and computes the group order `n` from `g`.
    pub fn new(g: Point<T>, ec: Box<dyn EllipticCurve<T>>) -> Self {
        // Check that the generator is valid.
        assert!(ec.is_valid(&g), "g must be a valid point on the elliptic curve");

        // Compute the group order by repeatedly adding g until the identity is reached.
        let n = ec.order(&g).expect("Unable to compute group order from g");

        DSA { g, n, ec }
    }

    pub fn gen_key(&self) -> (i32, Point<T>) {
        let mut rng_i = rand::rng();
        let n_i32 = self.n.to_i32().expect("n should be convertible to i32");
        let priv_gen = rng_i.random_range(1..n_i32);

        let point_pub = self.ec.mul(priv_gen, &self.g.clone());

        (priv_gen, point_pub)
    }
}
