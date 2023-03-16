use crate::field::GaloisField;
use crate::field_element::FieldElement;
use crate::utils;

use std::cmp;
use std::ops;

#[derive(Debug, Clone)]
pub struct Polynomial<'a> {
    coeffs: Vec<FieldElement<'a>>,
    field: &'a GaloisField,
    var: String,
}

impl<'a> Polynomial<'a> {
    pub fn x(field: &'a GaloisField) -> Polynomial<'a> {
        Polynomial {
            coeffs: vec![field.zero(), field.one()],
            field,
            var: "x".to_string(),
        }
    }

    pub fn new(coeffs: Vec<FieldElement<'a>>, field: &'a GaloisField, var: &str) -> Polynomial<'a> {
        let field_ptr = Polynomial::get_field_ptr(&coeffs);

        for el in coeffs.iter() {
            if el.field as *const _ != field_ptr {
                panic!(
                    "Not all coefficients of the constructing 
                polynomial are lay in the same field!"
                );
            }
        }

        let zero = field.zero();

        Polynomial {
            coeffs: utils::remove_trailing_elements(coeffs, zero),
            field,
            var: var.to_string(),
        }
    }

    pub fn get_field_ptr(coeffs: &Vec<FieldElement<'a>>) -> *const GaloisField {
        match coeffs.first() {
            Some(element) => element.field as *const _,
            _ => std::ptr::null() as *const _,
        }
    }
}

impl<'a> cmp::PartialEq<Polynomial<'a>> for Polynomial<'a> {
    fn eq(&self, other: &Polynomial<'a>) -> bool {
        self.coeffs == other.coeffs
    }
}

impl<'a> ops::Add<Polynomial<'a>> for Polynomial<'a> {
    type Output = Polynomial<'a>;

    fn add(self, rhs: Polynomial<'a>) -> Self::Output {
        if self.field as *const _ != rhs.field as *const _ {
            panic!("Polynomials are biult over different fields!");
        }

        if self.var != rhs.var {
            panic!("Polynomials have different variable names!");
        }

        Polynomial {
            coeffs: utils::zip_longest_with_op(
                &self.coeffs,
                &rhs.coeffs,
                |a, b| a + b,
                self.field.zero(),
            ),
            field: self.field,
            var: self.var,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::slice::RChunks;

    use super::Polynomial;
    use crate::field::GaloisField;

    #[test]
    fn init_test() {
        let field = GaloisField::predef();

        assert_eq!(
            Polynomial::x(&field),
            Polynomial {
                coeffs: vec![field.zero(), field.one()],
                field: &field,
                var: "x".to_string()
            }
        );
    }

    #[test]
    fn add_test() {
        let field = GaloisField::predef();

        let empty_poly = Polynomial::new(vec![], &field, "x");
        let some_poly = Polynomial::new(
            vec![
                field.new_element(1),
                field.new_element(2),
                field.new_element(3),
            ],
            &field,
            "x",
        );
        let other_poly = Polynomial::new(
            vec![
                field.new_element(2),
                field.new_element(3),
                field.new_element(4),
                field.new_element(5),
                field.new_element(6),
            ],
            &field,
            "x",
        );

        let sum1 = Polynomial::new(
            vec![
                field.new_element(1),
                field.new_element(2),
                field.new_element(3),
            ],
            &field,
            "x",
        );

        let sum2 = Polynomial::new(
            vec![
                field.new_element(3),
                field.new_element(5),
                field.new_element(7),
                field.new_element(5),
                field.new_element(6),
            ],
            &field,
            "x",
        );

        assert_eq!(empty_poly.clone() + empty_poly.clone(), empty_poly);
        assert_eq!(empty_poly.clone() + some_poly.clone(), sum1);
        assert_eq!(some_poly.clone() + empty_poly, sum1);
        assert_eq!(some_poly.clone() + other_poly.clone(), sum2);
        assert_eq!(other_poly + some_poly, sum2);
    }
}
