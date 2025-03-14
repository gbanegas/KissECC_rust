use KissECC::ecc::{EllipticCurve, Point};
use KissECC::weierstrass_ecc::WeierstrassECC;



fn main() {

        // For demonstration, using i32 as T (note: in cryptographic code, T will usually be a big integer type)
        let ecc = WeierstrassECC { a: 2, b: 3, q: 17 };
        let p = Point { x: 5, y: 1 , z: 0};
        println!("Is point valid? {}", ecc.is_valid(&p));

}
