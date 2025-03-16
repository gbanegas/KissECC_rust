use num_traits::{Zero, One, FromPrimitive, ToPrimitive};
use num_integer::Integer;
use std::ops::{Add, Sub, Mul, Rem, Div, BitAnd, Shr};
use crate::ecc::{EllipticCurve,};
use crate::point::Point;
use crate::utils::Utils;


/// Edwards curve represented by the equation:
///     a*x² + y² = 1 + d*x²*y²   (mod q)
/// with identity element (0, 1).
pub struct EdwardsCurve<T> {
    pub a: T,
    pub d: T,
    pub q: T,
    /// i = 2^((q-1)/4) mod q used for x recovery.
    pub i: T,
    /// The identity (zero) point: (0, 1).
    pub zero: Point<T>,
}

impl<T> EdwardsCurve<T>
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
    + Rem<Output = T>
    + Div<Output = T>
    + Copy,
{
    /// Creates a new EdwardsCurve given the parameters a, d and the prime modulus q.
    /// Also computes i = 2^((q-1)/4) mod q and sets the identity point to (0, 1).
    pub fn new(a: T, d: T, q: T) -> Self {
        let one = T::one();
        let two = T::from(2u8);
        let four = T::from(4u8);
        let q_minus_one = q.clone() - one.clone();
        let exp = ((q_minus_one) / four)
            .to_u32()
            .expect("Failed to convert exponent to u32");
        let i = Utils::modpow(two, exp, q.clone());
        let zero = Point {
            x: T::zero(),
            y: one.clone(),
            z: T::zero(),
        };
        EdwardsCurve { a, d, q, i: i, zero }
    }

    /// Given a y-coordinate, recover the corresponding x-coordinate.
    ///
    /// The procedure is as follows:
    /// 1. Compute xx = (y² − 1) / (d*y² + 1) mod q.
    /// 2. Compute x = xx^((q+3)/8) mod q.
    /// 3. If x² ≠ xx mod q then set x = x*i mod q.
    /// 4. Finally, ensure that x is “even” (if not, replace x with q − x).
    pub fn xrecover(&self, y: T) -> Result<T, &'static str> {
        let one = T::one();
        let y2 = (y.clone() * y.clone()) % self.q.clone();
        let numerator = (y2.clone() - one.clone()) % self.q.clone();
        let denominator = (self.d.clone() * y2.clone() + one.clone()) % self.q.clone();
        let inv_denominator = Utils::mod_inv(denominator, self.q.clone())?;
        let xx = (numerator * inv_denominator) % self.q.clone();
        // Compute exponent (q+3)/8 as a u32.
        let exp = ((self.q.clone() + T::from(3u8)) / T::from(8u8))
            .to_u32()
            .expect("Failed to convert exponent");
        let mut x = Utils::modpow(xx.clone(), exp, self.q.clone());
        if ((x.clone() * x.clone()) - xx.clone()) % self.q.clone() != T::zero() {
            x = (x * self.i.clone()) % self.q.clone();
        }
        // Ensure x is “even”. Here, we check if x mod 2 is nonzero.
        if (x.clone() % T::from(2u8)) != T::zero() {
            x = self.q.clone() - x;
        }
        Ok(x)
    }
}

impl<T> EllipticCurve<T> for EdwardsCurve<T>
where
    T: Zero
    + std::fmt::Display
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
    + Shr<u32, Output = T> + std::fmt::Display,
{
    /// Checks if a given point (x, y) satisfies the Edwards curve equation:
    ///     a*x² + y² = 1 + d*x²*y² (mod q)
    fn is_valid(&self, point: &Point<T>) -> bool {
        let x = point.x.clone();
        let y = point.y.clone();
        let left = (self.a.clone() * x.clone() * x.clone() + y.clone() * y.clone()) % self.q.clone();
        let one = T::one();
        let right = (one + self.d.clone() * x.clone() * x * y.clone() * y) % self.q.clone();
        left == right
    }

    /// Given a y-coordinate, returns the two points on the curve with that y value.
    /// Uses `xrecover` to compute the corresponding x-coordinate.
    fn at(&self, y: T) -> Result<(Point<T>, Point<T>), &'static str> {
        let x = self.xrecover(y.clone())?;
        let neg_x = if x == T::zero() {
            x.clone()
        } else {
            self.q.clone() - x.clone()
        };
        Ok((
            Point { x: x.clone(), y: y.clone(), z: T::zero() },
            Point { x: neg_x, y, z: T::zero() },
        ))
    }

    /// Edwards addition.
    ///
    /// Given two points p = (x₁, y₁) and q = (x₂, y₂), their sum R = (x₃, y₃)
    /// is given by:
    ///   x₃ = (x₁*y₂ + x₂*y₁) / (1 + d*x₁*x₂*y₁*y₂)
    ///   y₃ = (y₁*y₂ + x₁*x₂) / (1 − d*x₁*x₂*y₁*y₂)
    /// where divisions are computed as multiplication by the modular inverse.
    fn add(&self, p: &Point<T>, q: &Point<T>) -> Point<T> {
        let one = T::one();
        let x1 = p.x.clone();
        let y1 = p.y.clone();
        let x2 = q.x.clone();
        let y2 = q.y.clone();
        let q = self.q.clone();

        let denom1 = (one.clone() + self.d.clone() * x1.clone() * x2.clone() * y1.clone() * y2.clone()) % q.clone();
        let denom2 = (one.clone() - self.d.clone() * x1.clone() * x2.clone() * y1.clone() * y2.clone()) % q.clone();
        let inv_denom1 = Utils::mod_inv(denom1, q.clone()).expect("Inverse exists in add");
        let inv_denom2 = Utils::mod_inv(denom2, q.clone()).expect("Inverse exists in add");

        let x3 = ((x1 * y2) + (x2 * y1)) * inv_denom1 % q.clone();
        let y3 = ((y1 * y2) + (x1 * x2)) * inv_denom2 % q.clone();
        Point { x: x3, y: y3, z: T::zero() }
    }

    /// Point doubling: simply adds the point to itself.
    fn double(&self, p: &Point<T>) -> Point<T> {
        self.add(p, p)
    }

    /// Scalar multiplication using the double-and-add algorithm.
    ///
    /// This method multiplies the point `p` by the scalar `n` in O(log n) steps.
    fn mul(&self, mut n: T, p: &Point<T>) -> Point<T> {
        let zero_point = self.zero.clone();
        let mut r = zero_point;
        let mut m2 = p.clone();
        let one = T::one();
        while n > T::zero() {
            if (n.clone() & one.clone()) == one.clone() {
                r = self.add(&r, &m2);
            }
            n = n >> 1;
            m2 = self.add(&m2, &m2);
        }
        r
    }

    /// Computes the order of a point g by repeatedly adding it until the identity is reached.
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

    /// Returns a string representation of the Edwards curve equation.
    ///
    /// For example: "(a*x² + y² = 1 + d*x²*y²) mod q"
    fn display(&self) -> String {
        format!("({}*x² + y² = 1 + {}*x²*y²) mod {}", self.a, self.d, self.q)
    }
}

