pub mod field;
pub mod field_element;
pub mod polynomial;
pub mod utils;

fn main() {
    println!("Hello, world!");

    let field = crate::galois_field!();

    let el1 = field.new_element(10);
    let el2 = field.new_element(12);

    println!("sum is {:?}", el1 + el2);
}
