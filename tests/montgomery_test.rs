#[cfg(test)]
mod tests {
    use KissECC::ecc::{EllipticCurve};
    use KissECC::point::Point;
    use KissECC::montgomery_curve::MontgomeryCurve;

    #[test]
    fn test_is_valid() {
        // Example parameters (A, B, q, order) for demonstration.
        let ecc = MontgomeryCurve::new(2, 3, 17, 19);
        let identity = Point { x: 0, y: 1, z: 0 };
        assert_eq!(ecc.is_valid(&identity), true);
        // Further validity tests would require known finite points on this curve.
    }

    #[test]
    fn test_add() {
        let ecc = MontgomeryCurve::new(2, 3, 17, 19);
        let identity = Point { x: 0, y: 1, z: 0 };
        // Choose a sample point (the following values are for illustration).
        let p = Point { x: 5, y: 8, z: 1 };
        let sum = ecc.add(&identity, &p);
        assert_eq!(sum, p);
    }

    #[test]
    fn test_mul() {
        let ecc = MontgomeryCurve::new(2, 3, 17, 19);
        let identity = Point { x: 0, y: 1, z: 0 };
        let p = Point { x: 5, y: 8, z: 1 };

        let r0 = ecc.mul(0, &p);
        assert_eq!(r0, identity);
    }
}
