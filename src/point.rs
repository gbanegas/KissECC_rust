use num_traits::{Zero, One};
use std::ops::{Rem};
use crate::utils::Utils;

/// A simple point structure on a curve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Point<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Point<T>
where
    T: Clone
    + Copy
    + One
    + Zero
    + PartialEq
    + Rem<Output = T>
    + From<u8>
    + std::cmp::PartialOrd
    + std::ops::Sub
    + std::ops::Sub<Output = T>
    + std::ops::Div<Output = T>
    + std::fmt::Debug,

{
    pub fn eq_affine(&self, other: &Self, q: T) -> bool {
        // If one point is the identity (conventionally z=0), check using the chosen identity representation.
        if self.z == T::zero() || other.z == T::zero() {
            // Compare based on your identity representation.
            return self == other;
        }
        // Otherwise, normalize (i.e., compare x/z and y/z).
        let inv_self_z = Utils::mod_inv(self.z, q.clone()).unwrap();
        let inv_other_z = Utils::mod_inv(other.z, q.clone()).unwrap();
        let x1 = (self.x.clone() * inv_self_z) % q.clone();
        let y1 = (self.y.clone() * inv_self_z) % q.clone();
        let x2 = (other.x.clone() * inv_other_z) % q.clone();
        let y2 = (other.y.clone() * inv_other_z) % q;
        x1 == x2 && y1 == y2
    }
}

