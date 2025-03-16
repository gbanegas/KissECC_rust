
use num_traits::{Zero, One, FromPrimitive, ToPrimitive};
use std::ops::{Mul, Rem, Add, Sub, BitAnd, Shr, Div};
use std::cmp::PartialEq;
use num_integer::Integer;
use crate::ecc::{EllipticCurve};
use crate::point::Point;
use crate::utils::Utils;

pub struct WeierstrassECC<T> {
    //elliptic curve as: (y**2 = x**3 + a * x + b) mod q
    pub a: T,
    pub b: T,
    pub q: T,
    // You can include additional parameters as needed.
}


impl<T> WeierstrassECC<T>
where
    T: One
    + From<u8>
    + Clone
    + PartialEq
    + PartialOrd
    + FromPrimitive
    + ToPrimitive
    + Integer
    + Add<Output = T>
    + Sub<Output = T>
    + Mul<Output = T>
    + Rem<Output = T>
    + Div<Output = T>
    + Copy,
{
    /// Creates a new Weierstrass Curve
    ///

    pub fn new(a: T, b: T, q: T) -> Self {
        assert!(q > T::from(2u8));
        assert!(a != T::zero());
        assert!(b != T::zero());
        assert!(a != b);
        WeierstrassECC { a, b, q }
    }
}

impl<T> EllipticCurve<T> for WeierstrassECC<T>
where
    T: Zero
    + From<u8>
    + One
    + Clone
    + PartialEq
    + PartialOrd
    + FromPrimitive
    + ToPrimitive
    + Integer
    + Add<Output = T>
    + Sub<Output = T>
    + Mul<Output = T>
    + Copy
    + Rem<Output = T> + std::fmt::Display
    + BitAnd<Output = T>  // Required for n & 1
    + Shr<u32, Output = T>, // Required for n >> 1,
{

    fn is_valid(&self, p: &Point<T>) -> bool {
        // Check if p is the "zero" point.
        if p.x == T::zero() && p.y == T::zero() {
            return true;
        }
        // For example, for a curve defined by y^2 = x^3 + ax + b (mod q):
        let lhs = (p.y.clone() * p.y.clone()) % self.q.clone();
        let rhs = ((p.x.clone() * p.x.clone() * p.x.clone())
            + (self.a.clone() * p.x.clone())
            + self.b.clone())
            % self.q.clone();
        lhs == rhs
    }

    fn at(&self, x: T) -> Result<(Point<T>, Point<T>), &'static str> {
        // Ensure x < q; in practice, you'd want to handle this with proper error management.
        if x >= self.q.clone() {
            return Err("x must be less than q");
        }
        // Compute y^2 = x^3 + ax + b mod q.
        let rhs = ((x.clone() * x.clone() * x.clone())
            + (self.a.clone() * x.clone())
            + self.b.clone())
            % self.q.clone();

        // Use your modular square root function here.
        // For demonstration, we'll assume a function `mod_sqrt` exists:
        match Utils::tonelli_shanks(rhs, self.q.clone()) {
            Some(y) => {
                // Assuming modular negation is simply q - y
                let neg_y = (self.q.clone() - y.clone()) % self.q.clone();
                Ok((Point { x: x.clone(), y, z: T::zero() }, Point { x, y: neg_y, z: T::zero()}))
            },
            None => Err("No square root found"),
        }
    }

    fn add(&self, _p: &Point<T>, _q: &Point<T>) -> Point<T> {
        assert_eq!(self.is_valid(_p), true);
        assert_eq!(self.is_valid(_q), true);


        if _p.x == T::zero() && _p.y == T::zero() && _p.z == T::zero() {
            return _q.clone();
        }
        if _q.x == T::zero() && _q.y == T::zero() && _q.z == T::zero() {
            return _p.clone();
        }
        let mut res = Point { x: T::zero(), y: T::zero() , z: T::zero()};
        if _p.x == _q.x && (_p.y != _q.y || _p.y == T::zero()){
            return  res.clone();
        }
        let l;
        if _p.x == _q.x{
            // p1 + p1: use tangent line of p1 as (p1,p1) line

            let  two = T::from(2u8);
            let inv_val = Utils::mod_inv(two * _p.y.clone(), self.q.clone())
                .expect("Inverse should exist");
            let three = T::from(3u8);
            l = (((three.clone() * _p.x.clone() * _q.x.clone()) + self.a.clone()) * inv_val) % self.q.clone();
        }
        else {
            let tmp = (_q.x.clone() - _p.x.clone())%self.q.clone();
            let inv = Utils::mod_inv(tmp, self.q.clone()).expect("Inverse should exist");
            l = (_q.y - _p.y.clone())*inv.clone()%self.q.clone();
        }
        res.x = ((l*l)-_p.x - _q.x.clone()) % self.q.clone();
        res.y = (l*(_p.x-res.x) - _q.y.clone()) % self.q.clone();

        res
    }

    fn double(&self, _p: &Point<T>) -> Point<T> {
        // Implementation of point doubling goes here.
        unimplemented!()
    }

    fn mul(&self,  n: T, p: &Point<T>) -> Point<T> {
        // Define the zero (identity) point.
        let zero_point = Point {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        };

        let mut n_c = n.clone();

        // Start with r as the identity and m2 as p.
        let mut r = zero_point.clone();
        let mut m2 = p.clone();

        // Loop until n is reduced to zero.
        while n_c > T::zero() {
            // If the least significant bit of n is 1, add m2 into the result.
            if (n_c & T::one()) == T::one() {
                r = self.add(&r, &m2);
            }
            // Right shift n by 1 (divide by 2).
            n_c = n_c >> 1;
            // Double the point m2.
            m2 = self.add(&m2, &m2);
        }
        r.x = ((r.x % self.q.clone()) + self.q.clone()) % self.q.clone();
        r.y = ((r.y % self.q.clone()) + self.q.clone()) % self.q.clone();

        r
    }



    fn order(&self, g: &Point<T>) -> Result<T, &'static str> {
        // Define the identity (zero) point.
        let zero = Point {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        };

        // Ensure g is valid and not the zero point.
        if !self.is_valid(g) || *g == zero {
            return Err("Invalid point");
        }

        let mut order = T::one();
        let mut current = g.clone();
        // Repeatedly add g until we get the zero point.
        while current != zero {
            order = order + T::one();
            current = self.add(&current, g);
            // If we exceed the modulus, something is wrong.
            if order > self.q {
                return Err("Point order not found within group bounds");
            }
        }
        Ok(order)
    }

    fn display(&self) -> String {
        format!("(y**2 = x**3 + {} * x + {}) mod {}", self.a, self.b, self.q)
    }
}
