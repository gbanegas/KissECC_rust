#[cfg(test)]
mod tests {
    use KissECC::ecc::{EllipticCurve};
    use KissECC::point::Point;
    use KissECC::weierstrass_ecc::WeierstrassECC;

    #[test]
    fn test_is_valid() {
        let ecc = WeierstrassECC::new(  2, 3, 17 );
        // A point that is not on the curve.
        let invalid = Point { x: 5, y: 1, z: 0 };
        assert_eq!(ecc.is_valid(&invalid), false);
        // A known valid point on y² = x³ + 2x + 3 mod 17.
        let valid = Point { x: 5, y: 6, z: 0 };
        assert_eq!(ecc.is_valid(&valid), true);
        // The "zero" point (identity) is defined as (0, 0, 0) in this implementation.
        let identity = Point { x: 0, y: 0, z: 0 };
        assert_eq!(ecc.is_valid(&identity), true);
    }

    #[test]
    fn test_add() {
        let ecc = WeierstrassECC::new(  2, 3, 17 );
        let identity = Point { x: 0, y: 0, z: 0 };
        let p = Point { x: 5, y: 6, z: 0 };
        // Adding the identity should return the same point.
        let sum = ecc.add(&p, &identity);
        assert_eq!(sum, p);
    }

    #[test]
    fn test_mul() {
        let ecc = WeierstrassECC::new(  2, 3, 17 );
        let identity = Point { x: 0, y: 0, z: 0 };
        let p = Point { x: 5, y: 6, z: 1 };
        let p2 = Point { x: 15, y: 12, z: 1 };

        // Multiplying by 0 should yield the identity.
        let r0 = ecc.mul(0, &p);
        assert_eq!(r0, identity);

        // Multiplying by 1 should yield the same point.
        let r1 = ecc.mul(1, &p);
        assert_eq!(r1, p);

        let r2 = ecc.mul(2, &p);

        assert_eq!( r2.eq_affine(&p2, ecc.q), true);
    }
}
