
#[allow(dead_code)]


#[cfg(test)]
mod tests {
    use KissECC::ecc::{EllipticCurve, Point};
    use KissECC::twisted_curve::TwistedCurve;

    #[test]
    fn test_is_at(){

        let ecc = TwistedCurve::new( 2,  3,  17, 0 );
        let p = Point { x: 5, y: 1 , z: 0};
        assert_eq!(ecc.is_valid(&p), false);
        let p2 = Point { x: 1, y: 3 , z: 0};
        assert_eq!(ecc.is_valid(&p2), true);
    }

    fn test_add(){
    }

    fn test_mul(){
    }


}
