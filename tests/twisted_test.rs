#[cfg(test)]
mod tests {
    use KissECC::ecc::{EllipticCurve};
    use KissECC::point::Point;
    use KissECC::twisted_curve::TwistedCurve;

    #[test]
    fn test_is_valid() {
        // Example parameters: a, b, q, order.
        let ecc = TwistedCurve::new(2, 3, 17, 19);
        let invalid = Point { x: 5, y: 1, z: 0 };
        assert_eq!(ecc.is_valid(&invalid), false);
        let valid = Point { x: 1, y: 3, z: 0 };
        assert_eq!(ecc.is_valid(&valid), true);
    }

    #[test]
    fn test_add() {
        let ecc = TwistedCurve::new(2, 3, 17, 19);
        let identity = ecc.zero.clone();
        let p = Point { x: 1, y: 3, z: 0 };
        let sum = ecc.add(&identity, &p);
        assert_eq!(sum, p);
    }

    #[test]
    fn test_mul() {
        let ecc = TwistedCurve::new(2, 3, 17, 19);
        let identity = ecc.zero.clone();
        let p = Point { x: 1, y: 3, z: 0 };

        let r0 = ecc.mul(0, &p);
        assert_eq!(r0, identity);

        let r1 = ecc.mul(1, &p);
        assert_eq!(r1, p);
    }
}
