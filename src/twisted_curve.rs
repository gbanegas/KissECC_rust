use num_traits::{Zero, One, FromPrimitive, ToPrimitive};
use num_integer::Integer;
use std::ops::{Add, Sub, Mul, Rem, Div, BitAnd, Shr};
use crate::ecc::{EllipticCurve};
use crate::point::Point;
use crate::utils::Utils;

/// Twisted Edwards curve defined by the equation:
///     a*x² + y² = 1 + b*x²*y²   (mod q)
///
/// The curve stores:
/// - a, b: parameters (with a ≠ 0, b ≠ 0 and a ≠ b),
/// - q: a prime number > 2,
/// - I: computed as 2^((q-1)/4) mod q (used in x recovery),
/// - zero: the identity element, here (0, 1),
/// - order: the group order.

#[allow(non_snake_case)]
pub struct TwistedCurve<T> {
    pub a: T,
    pub b: T,
    pub q: T,
    pub I: T,
    pub zero: Point<T>,
    pub order: T,
}

impl<T> TwistedCurve<T>
where
    T: Zero
    + One
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
    /// Creates a new twisted Edwards curve.
    ///
    /// It asserts that q > 2, a and b are nonzero and different, and that an order is provided.
    /// Also computes I = 2^((q-1)/4) mod q and sets the identity (zero) point as (0, 1).
    pub fn new(a: T, b: T, q: T, order: T) -> Self {
        assert!(q > T::from(2u8));
        assert!(a != T::zero());
        assert!(b != T::zero());
        assert!(a != b);
        let one = T::one();
        let two = T::from(2u8);
        // Compute exponent = (q-1)/4 as a u32.
        let exp = ((q.clone() - one.clone()) / T::from(4u8))
            .to_u32()
            .expect("Failed to convert exponent");
        let i = Utils::modpow(two, exp, q.clone());
        let zero = Point {
            x: T::zero(),
            y: one.clone(),
            z: T::zero(),
        };
        TwistedCurve { a, b, q, I: i, zero, order }
    }

    /// Recovers the x-coordinate corresponding to a given y-coordinate.
    ///
    /// The procedure is as follows:
    /// 1. Compute: xx = (y² − 1) / (b*y² + 1)  (using a modular inverse).
    /// 2. Compute: x = xx^((q+3)/8) mod q.
    /// 3. If \(x^2 \not\equiv xx \bmod q\) then adjust \(x = x * I \bmod q\).
    /// 4. Finally, ensure that \(x\) is “even” (if not, replace \(x\) with \(q - x\)).
    pub fn xrecover(&self, y: T) -> Result<T, &'static str> {
        let one = T::one();
        let y2 = (y.clone() * y.clone()) % self.q.clone();
        let numerator = (y2.clone() - one.clone()) % self.q.clone();
        let denominator = (self.b.clone() * y2.clone() + one.clone()) % self.q.clone();
        let inv_den = Utils::mod_inv(denominator, self.q.clone())?;
        let xx = (numerator * inv_den) % self.q.clone();
        let exp = ((self.q.clone() + T::from(3u8)) / T::from(8u8))
            .to_u32()
            .expect("Failed to convert exponent");
        let mut x = Utils::modpow(xx.clone(), exp, self.q.clone());
        if ((x.clone() * x.clone()) - xx.clone()) % self.q.clone() != T::zero() {
            x = (x * self.I.clone()) % self.q.clone();
        }
        if (x.clone() % T::from(2u8)) != T::zero() {
            x = self.q.clone() - x;
        }
        Ok(x)
    }

    /// A version of Edwards addition tailored for the twisted curve.
    ///
    /// Given two points \(p = (x_1, y_1)\) and \(q = (x_2, y_2)\), their sum \(R = (x_3, y_3)\)
    /// is defined by:
    ///
    /// \[
    /// x_3 = \frac{x_1y_2 + x_2y_1}{1 + b x_1x_2y_1y_2} \quad\text{and}\quad
    /// y_3 = \frac{y_1y_2 - a x_1x_2}{1 - b x_1x_2y_1y_2} \quad \text{(mod }q\text{)}
    /// \]
    ///
    /// Division is performed by multiplying by the modular inverse.
    pub fn edwards_add(&self, p: &Point<T>, _q: &Point<T>) -> Point<T> {
        let one = T::one();
        let q = self.q.clone();
        let x1 = p.x.clone();
        let y1 = p.y.clone();
        let x2 = _q.x.clone();
        let y2 = _q.y.clone();
        let factor = (self.b.clone() * x1.clone() * x2.clone() * y1.clone() * y2.clone()) % q.clone();
        let denom1 = (one.clone() + factor.clone()) % q.clone();
        let denom2 = (one.clone() - factor) % q.clone();
        let inv_denom1 = Utils::mod_inv(denom1, q.clone()).expect("Inverse exists in edwards_add");
        let inv_denom2 = Utils::mod_inv(denom2, q.clone()).expect("Inverse exists in edwards_add");
        let x3 = (((x1 * y2.clone()) + (x2 * y1.clone())) * inv_denom1) % q.clone();
        let y3 = (((y1 * y2) - (self.a.clone() * x1 * x2)) * inv_denom2) % q.clone();
        Point { x: x3, y: y3, z: T::zero() }
    }
}

impl<T> EllipticCurve<T> for TwistedCurve<T>
where
    T: Zero
    + One
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
    + Copy
    + BitAnd<Output = T>
    + std::fmt::Display
    + Shr<u32, Output = T>,
{
    /// Checks whether a point \(P = (x, y)\) satisfies the twisted Edwards curve equation:
    ///     a*x² + y² = 1 + b*x²*y²  (mod q)
    ///
    /// (Note: the Python version checks that \(-x^2 + y^2 - 1 - b*x^2*y^2 \equiv 0\), which is equivalent when \(a=-1\).)
    fn is_valid(&self, point: &Point<T>) -> bool {
        let x = point.x.clone();
        let y = point.y.clone();
        let left = (T::zero() - (x.clone() * x.clone()) + (y.clone() * y.clone()) - T::one()
            - (self.b.clone() * x.clone() * x.clone() * y.clone() * y.clone()))
            % self.q.clone();
        left == T::zero()
    }

    /// Given a y-coordinate, returns the two corresponding points on the curve by recovering \(x\).
    fn at(&self, y: T) -> Result<(Point<T>, Point<T>), &'static str> {
        let x = self.xrecover(y.clone())?;
        let neg_x = if x == T::zero() { x.clone() } else { self.q.clone() - x.clone() };
        Ok((
            Point { x: x.clone(), y: y.clone(), z: T::zero() },
            Point { x: neg_x, y, z: T::zero() },
        ))
    }

    /// Adds two points on the curve.
    ///
    /// If either point is the identity, returns the other.
    /// Otherwise, uses the specialized Edwards addition for the twisted curve.
    fn add(&self, p: &Point<T>, q: &Point<T>) -> Point<T> {
        if *p == self.zero {
            return q.clone();
        }
        if *q == self.zero {
            return p.clone();
        }
        self.edwards_add(p, q)
    }

    /// Doubles a point by adding it to itself.
    fn double(&self, p: &Point<T>) -> Point<T> {
        self.add(p, p)
    }

    /// Scalar multiplication using the double-and-add algorithm.
    fn mul(&self, mut n: i32, p: &Point<T>) -> Point<T> {
        let zero_point = self.zero.clone();
        let mut r = zero_point;
        let mut m2 = p.clone();
        while n >0 {
            if (n.clone() & 1) == 1 {
                r = self.add(&r, &m2);
            }
            n = n >> 1;
            m2 = self.add(&m2, &m2);
        }
        r
    }

    /// Computes the order of a point by repeatedly adding it until the identity is reached.
    fn order(&self, g: &Point<T>) -> Result<T, &'static str> {
        let mut order = T::one();
        let mut current = g.clone();
        let zero_point = self.zero.clone();
        while current != zero_point {
            order = order + T::one();
            current = self.add(&current, g);
            if order > self.q {
                return Err("Order not found within group bounds");
            }
        }
        Ok(order)
    }

    /// Returns a string representation of the twisted Edwards curve.
    fn display(&self) -> String {
        format!("({}*x² + y² = 1 + {}*x²*y²) mod {}", self.a, self.b, self.q)
    }
}
