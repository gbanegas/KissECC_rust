use num_traits::{Zero, One, FromPrimitive, ToPrimitive};
use num_integer::Integer;
use std::ops::{Add, Sub, Mul, Rem, Div, BitAnd, Shr};
use crate::ecc::{EllipticCurve};
use crate::point::Point;
use crate::utils::Utils;

/// Montgomery curve defined by:
///      B * y^2 = x^3 + A * x^2 + x   (mod q)
///
/// We include the parameters A, B, the prime modulus q, and the group order.
/// The identity (zero) point is represented as (0, 1, 0).
///
#[allow(non_snake_case)]
pub struct MontgomeryCurve<T> {
    pub A: T,
    pub B: T,
    pub q: T,
    pub order: T,
    pub zero: Point<T>,
}

impl<T> MontgomeryCurve<T>
where
    T: Zero
    + One
    + std::fmt::Display
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
    /// Creates a new Montgomery curve with given parameters.
    /// Expects q to be a prime number > 2, and that the curve order is provided.
    #[allow(non_snake_case)]
    pub fn new(A: T, B: T, q: T, order: T) -> Self {
        assert!(q > T::from(2u8));
        // We assume A and B are nonzero.
        assert!(A != T::zero() && B != T::zero());
        // Define the identity element.
        let zero = Point {
            x: T::zero(),
            y: T::one(),
            z: T::zero(),
        };
        MontgomeryCurve { A, B, q, order, zero }
    }

    /// Helper: checks whether a given point is the identity.
    fn is_zero(&self, p: &Point<T>) -> bool {
        *p == self.zero
    }
}

impl<T> EllipticCurve<T> for MontgomeryCurve<T>
where
    T: Zero
    + One
    + Clone
    + From<u8>
    + std::fmt::Display
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
    + Shr<u32, Output = T>,
{
    /// A point is valid if it is the identity, or if it satisfies
    ///      B*y^2 = x^3 + A*x^2 + x  (mod q).
    fn is_valid(&self, p: &Point<T>) -> bool {
        if self.is_zero(p) {
            return true;
        }
        let x = p.x.clone();
        let y = p.y.clone();
        let left = (self.B * y.clone() * y.clone()) % self.q.clone();
        let right = (x.clone() * x.clone() * x.clone()
            + self.A * x.clone() * x.clone()
            + x) % self.q.clone();
        left == right
    }

    /// Given an x-coordinate, we cannot in general recover y uniquely on a Montgomery curve.
    /// Here we return an error.
    fn at(&self, _x: T) -> Result<(Point<T>, Point<T>), &'static str> {
        Err("Method 'at' is not implemented for Montgomery curves")
    }

    /// Adds two points p and Q.
    ///
    /// - If either point is the identity, returns the other.
    /// - If p and Q have the same x-coordinate but different y, returns the identity.
    /// - Otherwise, uses the chord-and-tangent formulas:
    ///
    /// For p != Q:
    ///   λ = (y₂ − y₁)/(x₂ − x₁)
    ///   x₃ = B*λ² − A − x₁ − x₂
    ///   y₃ = λ*(x₁ − x₃) − y₁
    ///
    /// For doubling (p == Q):
    ///   λ = (3*x₁² + 2*A*x₁ + 1)/(2*B*y₁)
    ///   x₃ = B*λ² − A − 2*x₁
    ///   y₃ = λ*(x₁ − x₃) − y₁
    fn add(&self, p: &Point<T>, _q: &Point<T>) -> Point<T> {
        // Handle identity cases.
        if self.is_zero(p) {
            return _q.clone();
        }
        if self.is_zero(_q) {
            return p.clone();
        }
        // If x1 == x2 and y1 != y2, then p + _q = 0.
        if p.x == _q.x && p.y != _q.y {
            return self.zero.clone();
        }
        // Determine whether we are doubling or adding distinct points.
        let (_lambda, x3, y3) = if p == _q {
            // Doubling: p = _q.
            // λ = (3*x₁² + 2*A*x₁ + 1)/(2*B*y₁)
            let numerator = (T::from(3u8) * p.x.clone() * p.x.clone()
                + (T::from(2u8) * self.A * p.x.clone())
                + T::one()) % self.q.clone();
            let denominator = (T::from(2u8) * self.B * p.y.clone()) % self.q.clone();
            let inv_den = Utils::mod_inv(denominator, self.q.clone())
                .expect("Denom invertible in doubling");
            let lambda = (numerator * inv_den) % self.q.clone();
            // x₃ = B*λ² − A − 2*x₁
            let x3 = (self.B * lambda.clone() * lambda.clone()
                - self.A - (T::from(2u8) * p.x.clone())) % self.q.clone();
            // y₃ = λ*(x₁ − x₃) − y₁
            let y3 = (lambda * (p.x.clone() - x3.clone()) - p.y.clone()) % self.q.clone();
            (lambda, x3, y3)
        } else {
            // Addition: p != _q.
            // λ = (y₂ − y₁)/(x₂ − x₁)
            let numerator = (_q.y.clone() - p.y.clone()) % self.q.clone();
            let denominator = (_q.x.clone() - p.x.clone()) % self.q.clone();
            let inv_den = Utils::mod_inv(denominator, self.q.clone())
                .expect("Denom invertible in addition");
            let lambda = (numerator * inv_den) % self.q.clone();
            // x₃ = B*λ² − A − x₁ − x₂
            let x3 = (self.B * lambda.clone() * lambda.clone()
                - self.A - p.x.clone() - _q.x.clone()) % self.q.clone();
            // y₃ = λ*(x₁ − x₃) − y₁
            let y3 = (lambda * (p.x.clone() - x3.clone()) - p.y.clone()) % self.q.clone();
            (lambda, x3, y3)
        };
        // Normalize (make sure the result is positive modulo q).
        Point {
            x: (x3 + self.q.clone()) % self.q.clone(),
            y: (y3 + self.q.clone()) % self.q.clone(),
            z: T::one(), // For finite points we set z = 1.
        }
    }

    /// Doubles a point using the formula described in `add`.
    fn double(&self, p: &Point<T>) -> Point<T> {
        self.add(p, p)
    }

    /// Scalar multiplication via the double-and-add algorithm.
    fn mul(&self, mut n: T, p: &Point<T>) -> Point<T> {
        let mut r = self.zero.clone();
        let mut m2 = p.clone();
        let one = T::one();
        while n > T::zero() {
            if (n.clone() & one.clone()) == one {
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
        while current != self.zero {
            order = order + T::one();
            current = self.add(&current, g);
            if order > self.q {
                return Err("Order not found within group bounds");
            }
        }
        Ok(order)
    }

    /// Returns a string representation of the Montgomery curve.
    fn display(&self) -> String {
        format!("({}*y² = x³ + {}*x² + x) mod {}", self.B, self.A, self.q)
    }
}
