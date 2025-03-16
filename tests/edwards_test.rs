// This is a really bad adding function, its purpose is to fail in this
// example.
#[allow(dead_code)]


#[cfg(test)]
mod tests {
    use KissECC::ecc::{EllipticCurve, Point};
    use KissECC::edwards_curve::EdwardsCurve;

    #[test]
    fn test_is_at(){

        let ecc = EdwardsCurve::new( 2,  3,  17 );
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
