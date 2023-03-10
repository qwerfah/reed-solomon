pub mod field;

use crate::field::field::*;

fn main() {
    println!("Hello, world!");

    let field = GaloisField::predef();

    let el1 = field.new_element(10);
    let el2 = field.new_element(12);

    println!("sum is {:?}", el1 + el2);
}
