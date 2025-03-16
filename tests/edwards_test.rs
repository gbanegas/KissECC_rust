#[cfg(test)]
mod tests {
    use KissECC::ecc::{EllipticCurve, Point};
    use KissECC::edwards_curve::EdwardsCurve;

    #[test]
    fn test_is_valid() {
        let ecc = EdwardsCurve::new(2, 3, 17);
        let invalid = Point { x: 5, y: 1, z: 0 };
        assert_eq!(ecc.is_valid(&invalid), false);
        let valid = Point { x: 1, y: 3, z: 0 };
        assert_eq!(ecc.is_valid(&valid), true);
    }

    #[test]
    fn test_add() {
        let ecc = EdwardsCurve::new(2, 3, 17);
        // The identity for Edwards is (0, 1, 0)
        let identity = Point { x: 0, y: 1, z: 0 };
        let p = Point { x: 1, y: 3, z: 0 };
        let sum = ecc.add(&identity, &p);
        assert_eq!(sum, p);
    }

    #[test]
    fn test_mul() {
        let ecc = EdwardsCurve::new(2, 3, 17);
        let identity = Point { x: 0, y: 1, z: 0 };
        let p = Point { x: 1, y: 3, z: 0 };

        let r0 = ecc.mul(0, &p);
        assert_eq!(r0, identity);

        let r1 = ecc.mul(1, &p);
        assert_eq!(r1, p);
    }
}
