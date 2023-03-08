use crate::field::field_element::FieldElement;

/// Finite field with order of 2^power. Cause we always works with a chunks of data
/// multiple of octet, the number of field elements always will be a power of 2.
#[derive(Debug)]
pub struct GaloisField {
    pub k_modulus: u32,
    pub generator_val: u32,

    pub zero: u32,
    pub one: u32,
}

impl GaloisField {
    pub fn predef() -> GaloisField {
        GaloisField {
            k_modulus: 3 * u32::pow(2, 30) + 1,
            generator_val: 5,
            zero: 0,
            one: 1,
        }
    }

    pub fn new_element<'a>(&'a self, element_val: u32) -> FieldElement {
        FieldElement {
            val: element_val % self.k_modulus,
            field: self,
        }
    }

    pub fn zero<'a>(&'a self) -> FieldElement {
        FieldElement {
            val: self.zero,
            field: self,
        }
    }

    pub fn one<'a>(&'a self) -> FieldElement {
        FieldElement {
            val: self.one,
            field: self,
        }
    }

    pub fn generator<'a>(&'a self) -> FieldElement {
        FieldElement {
            val: self.generator_val,
            field: self,
        }
    }
}
