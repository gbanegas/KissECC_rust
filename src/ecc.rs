use num_traits::{Zero, One, FromPrimitive, ToPrimitive};
use std::ops::{Add, Sub, Mul, Rem};
use num_integer::Integer;

/// A simple point structure on a curve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Point<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

/// A trait that defines the common operations for an elliptic curve.
pub trait EllipticCurve<T>
where
    T: Zero
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
    + Rem<Output = T>,
{


    /// Checks whether a given point is valid on the curve.
    fn is_valid(&self, point: &Point<T>) -> bool;

    /// Given an x-coordinate (with x < q), returns the two points on the curve,
    /// i.e. (x, y) and (x, -y), if a square root exists.
    ///
    /// # Errors
    ///
    /// Returns an error if no valid square root is found.
    fn at(&self, x: T) -> Result<(Point<T>, Point<T>), &'static str>;

    /// Adds two points on the curve.
    fn add(&self, p: &Point<T>, q: &Point<T>) -> Point<T>;

    /// Doubles a point on the curve.
    fn double(&self, p: &Point<T>) -> Point<T>;

    /// Multiplies a point by a scalar k (i.e. repeated addition).
    fn mul(&self, k: T, p: &Point<T>) -> Point<T>;

    /// Returns the order of the curve (or the group order).
    fn order(&self) -> T;
}
