use num_traits::{One, Zero, FromPrimitive, ToPrimitive};
use num_integer::Integer;
use std::ops::{Add, Div, Mul, Rem, Sub};

pub struct Utils;

impl Utils {




    /// Provided modular exponentiation function.
    pub fn modpow<T>(base: T, exp: u32, modulus: T) -> T
    where
        T: Clone + PartialEq + Zero + One + Mul<Output = T> + Rem<Output = T> + ToPrimitive,
    {
        let mut result = T::one();
        let mut b = base % modulus.clone();
        let mut e = exp;
        while e > 0 {
            if e % 2 == 1 {
                result = (result * b.clone()) % modulus.clone();
            }
            e /= 2;
            b = (b.clone() * b.clone()) % modulus.clone();
        }
        result
    }

    /// Tonelli–Shanks algorithm: given a quadratic residue `a` modulo a prime `p`,
    /// finds an `x` such that \( x^2 \equiv a \) (mod \( p \)) (if one exists).
    ///
    /// Returns `Some(x)` if a square root exists, otherwise `None`.
    pub fn tonelli_shanks<T>(a: T, p: T) -> Option<T>
    where
        T: Clone
        + PartialEq
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Rem<Output = T>
        + FromPrimitive
        + ToPrimitive
        + Integer,
    {
        // If "a" is zero, the square root is zero.
        if a == T::zero() {
            return Some(T::zero());
        }

        // Convert p to u32. (This implementation assumes p is small enough.)
        let p_u32 = p.to_u32().expect("p should convert to u32");

        // Check the Legendre symbol: a^{(p-1)/2} mod p should equal 1 for a quadratic residue.
        let exp = (p_u32 - 1) / 2;
        if Utils::modpow(a.clone(), exp, p.clone()) != T::one() {
            return None; // a is not a quadratic residue modulo p.
        }

        // Write p - 1 as q * 2^s with q odd.
        let mut q = p_u32 - 1;
        let mut s = 0;
        while q % 2 == 0 {
            q /= 2;
            s += 1;
        }

        // Find a quadratic non-residue z (i.e. one for which the Legendre symbol is -1).
        let mut z = T::from_u32(2).unwrap();
        while Utils::modpow(z.clone(), exp, p.clone()) == T::one() {
            z = z + T::one();
        }

        // Set initial values:
        // c = z^q mod p,
        // x = a^{(q+1)/2} mod p,
        // t = a^q mod p.
        let mut c = Utils::modpow(z, q, p.clone());
        let mut x = Utils::modpow(a.clone(), (q + 1) / 2, p.clone());
        let mut t = Utils::modpow(a, q, p.clone());
        let mut m = s;

        // Main loop: adjust x until t ≡ 1 (mod p).
        while t != T::one() {
            // Find the smallest i (0 < i < m) such that t^(2^i) ≡ 1 mod p.
            let mut i = 1;
            while i < m && Utils::modpow(t.clone(), 1 << i, p.clone()) != T::one() {
                i += 1;
            }
            if i == m {
                return None; // This should not happen if a square root exists.
            }
            // Compute b = c^(2^(m-i-1)) mod p.
            let b = Utils::modpow(c.clone(), 1 << (m - i - 1), p.clone());
            x = (x * b.clone()) % p.clone();
            t = (t * b.clone() * b.clone()) % p.clone();
            c = (b.clone() * b.clone()) % p.clone();
            m = i;
        }
        Some(x)
    }

    /// Computes the modular inverse of `a` modulo `q` using the Extended Euclidean Algorithm.
    /// Returns an error if the inverse does not exist.
    pub fn mod_inv<T>(a: T, q: T) -> Result<T, &'static str>
    where
        T: Copy
        + PartialEq
        + PartialOrd
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Rem<Output = T>,
    {
        let zero = T::zero();
        let one = T::one();

        let mut t = zero;
        let mut new_t = one;
        let mut r = q;
        let mut new_r = a;

        while new_r != zero {
            let quotient = r / new_r;
            let temp_t = new_t;
            new_t = t - quotient * new_t;
            t = temp_t;

            let temp_r = new_r;
            new_r = r - quotient * new_r;
            r = temp_r;
        }

        if r > one {
            return Err("Inverse does not exist because a and q are not coprime");
        }
        if t < zero {
            t = t + q;
        }
        Ok(t)
    }


}