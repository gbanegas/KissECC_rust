use KissECC::ecc::{EllipticCurve, Point};
use KissECC::twisted_curve::TwistedCurve;



fn main() {

        // For demonstration, using i32 as T (note: in cryptographic code, T will usually be a big integer type)
        let ecc = TwistedCurve::new( 2,  3,  17, 0 );
        let p = Point { x: 1, y: 3 , z: 0};
        println!("Is point valid? {}", ecc.is_valid(&p));

        println!("{}",ecc.display());

}
