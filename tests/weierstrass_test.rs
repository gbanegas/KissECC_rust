// This is a really bad adding function, its purpose is to fail in this
// example.
#[allow(dead_code)]


#[cfg(test)]
mod tests {
    use KissECC::ecc::{EllipticCurve, Point};
    use KissECC::weierstrass_ecc::WeierstrassECC;

    #[test]
    fn test_is_at(){

        let ecc = WeierstrassECC { a: 2, b: 3, q: 17 };
        let p = Point { x: 5, y: 1 , z: 0};
        assert_eq!(ecc.is_valid(&p), false);
        let p2 = Point { x: 5, y: 6 , z: 0};
        assert_eq!(ecc.is_valid(&p2), true);
        let p3 = Point { x: 0, y: 0 , z: 0};
        assert_eq!(ecc.is_valid(&p3), true);
    }

    fn test_add(){
    }

    fn test_mul(){
    }


}
