use num_traits::{Zero, One, FromPrimitive, ToPrimitive};
use std::ops::{Mul, Rem, Add, Sub, BitAnd, Shr, Div};
use num_integer::Integer;
use crate::ecc::{EllipticCurve};
use crate::point::Point;
use crate::utils::Utils;

pub struct WeierstrassECC<T> {
    // elliptic curve: y² = x³ + a * x + b mod q
    pub a: T,
    pub b: T,
    pub q: T,
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
    /// Creates a new Weierstrass curve.
    pub fn new(a: T, b: T, q: T) -> Self {
        assert!(q > T::from(2u8));
        assert!(a != T::zero());
        assert!(b != T::zero());
        // You might want additional assertions depending on your use case.
        WeierstrassECC { a, b, q }
    }

    /// Normalize a point.
    ///
    /// If the point is not the identity, adjust x and y modulo q and set z to 1.
    pub fn normalize(&self, mut p: Point<T>) -> Point<T> {
        // Define the identity point as (0,0,0)
        let identity = Point {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        };
        if p == identity {
            p
        } else {
            p.x = ((p.x % self.q.clone()) + self.q.clone()) % self.q.clone();
            p.y = ((p.y % self.q.clone()) + self.q.clone()) % self.q.clone();
            p.z = T::one();
            p
        }
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
    + Rem<Output = T>
    + std::fmt::Display
    + BitAnd<Output = T>  // for n & 1
    + Shr<u32, Output = T>, // for n >> 1,
{

    fn is_valid(&self, p: &Point<T>) -> bool {
        // Treat (0,0) as the identity.
        if p.x == T::zero() && p.y == T::zero() {
            return true;
        }
        let lhs = (p.y.clone() * p.y.clone()) % self.q.clone();
        let rhs = ((p.x.clone() * p.x.clone() * p.x.clone())
            + (self.a.clone() * p.x.clone())
            + self.b.clone())
            % self.q.clone();
        lhs == rhs
    }

    fn at(&self, x: T) -> Result<(Point<T>, Point<T>), &'static str> {
        if x >= self.q.clone() {
            return Err("x must be less than q");
        }
        let rhs = ((x.clone() * x.clone() * x.clone())
            + (self.a.clone() * x.clone())
            + self.b.clone())
            % self.q.clone();
        match Utils::tonelli_shanks(rhs, self.q.clone()) {
            Some(y) => {
                let neg_y = (self.q.clone() - y.clone()) % self.q.clone();
                // Return normalized points.
                Ok((
                    self.normalize(Point { x: x.clone(), y: y, z: T::one() }),
                    self.normalize(Point { x, y: neg_y, z: T::one() }),
                ))
            },
            None => Err("No square root found"),
        }
    }

    fn add(&self, _p: &Point<T>, _q: &Point<T>) -> Point<T> {
        assert_eq!(self.is_valid(_p), true);
        assert_eq!(self.is_valid(_q), true);

        // Define the identity point.
        let identity = Point { x: T::zero(), y: T::zero(), z: T::zero() };
        if *_p == identity {
            return _q.clone();
        }
        if *_q == identity {
            return _p.clone();
        }
        // If x coordinates are equal and y differ (or y is zero), result is identity.
        if _p.x == _q.x && (_p.y != _q.y || _p.y == T::zero()) {
            return identity;
        }
        let l;
        if _p.x == _q.x {
            // Tangent case.
            let two = T::from(2u8);
            let inv_val = Utils::mod_inv(two * _p.y.clone(), self.q.clone())
                .expect("Inverse should exist");
            let three = T::from(3u8);
            l = (((three * _p.x.clone() * _q.x.clone()) + self.a.clone()) * inv_val) % self.q.clone();
        } else {
            // Chord case.
            let tmp = (_q.x.clone() - _p.x.clone()) % self.q.clone();
            let inv = Utils::mod_inv(tmp, self.q.clone()).expect("Inverse should exist");
            l = ((_q.y - _p.y.clone()) * inv) % self.q.clone();
        }
        let mut res = Point { x: T::zero(), y: T::zero(), z: T::zero() };
        res.x = ((l * l) - _p.x - _q.x.clone()) % self.q.clone();
        res.y = (l * (_p.x - res.x.clone()) - _p.y.clone()) % self.q.clone();
        // Normalize the result so that nonzero points have z = 1.
        self.normalize(res)
    }

    fn double(&self, p: &Point<T>) -> Point<T> {
        // Define the identity point.
        let identity = Point {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        };

        // If p is the identity, return it.
        if *p == identity {
            return p.clone();
        }
        // If y == 0, doubling yields the identity.
        if p.y == T::zero() {
            return identity;
        }

        // Calculate lambda = (3*x^2 + a) / (2*y) mod q.
        let two = T::from(2u8);
        let three = T::from(3u8);
        let x_sq = p.x.clone() * p.x.clone();
        let numerator = (three * x_sq + self.a.clone()) % self.q.clone();
        let denominator = (two * p.y.clone()) % self.q.clone();
        let inv_den = Utils::mod_inv(denominator, self.q.clone())
            .expect("Inverse should exist for denominator in doubling");
        let lambda = (numerator * inv_den) % self.q.clone();

        // x3 = lambda^2 - 2*x
        let lambda_sq = lambda.clone() * lambda.clone();
        let x3 = (lambda_sq - two * p.x.clone()) % self.q.clone();
        // y3 = lambda*(x - x3) - y
        let y3 = (lambda * (p.x.clone() - x3.clone()) - p.y.clone()) % self.q.clone();

        // Construct the result with z = 1 (after normalization).
        let result = Point {
            x: x3,
            y: y3,
            z: T::one(),
        };
        self.normalize(result)
    }


    fn mul(&self, n: T, p: &Point<T>) -> Point<T> {
        let zero_point = Point {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        };
        let mut n_c = n.clone();
        let mut r = zero_point.clone();
        let mut m2 = p.clone();
        while n_c > T::zero() {
            if (n_c & T::one()) == T::one() {
                r = self.add(&r, &m2);
            }
            n_c = n_c >> 1;
            m2 = self.add(&m2, &m2);
        }
        // Normalize the resulting point.
        self.normalize(r)
    }

    fn order(&self, g: &Point<T>) -> Result<T, &'static str> {
        let identity = Point { x: T::zero(), y: T::zero(), z: T::zero() };
        if !self.is_valid(g) || *g == identity {
            return Err("Invalid point");
        }
        let mut order = T::one();
        let mut current = g.clone();
        while current != identity {
            order = order + T::one();
            current = self.add(&current, g);
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
