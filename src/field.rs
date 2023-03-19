use crate::field_element::FieldElement;

/// Finite field with order of 2^power. Cause we always works with a chunks of data
/// multiple of octet, the number of field elements always will be a power of 2.
#[derive(Debug)]
pub struct GaloisField {
    pub k_modulus: u64,
    pub generator_val: u64,

    pub zero: u64,
    pub one: u64,
}

#[macro_export]
macro_rules! galois_field {
    () => {
        $crate::field::GaloisField {
            k_modulus: 3 * u64::pow(2, 30) + 1,
            generator_val: 5,
            zero: 0,
            one: 1,
        }
    };
}

impl GaloisField {
    pub fn new_element(&self, element_val: i128) -> FieldElement {
        FieldElement {
            val: element_val.rem_euclid(self.k_modulus as i128) as u64,
            field: self,
        }
    }

    pub fn zero(&'_ self) -> FieldElement {
        FieldElement {
            val: self.zero,
            field: self,
        }
    }

    pub fn one(&'_ self) -> FieldElement {
        FieldElement {
            val: self.one,
            field: self,
        }
    }

    pub fn generator(&'_ self) -> FieldElement {
        FieldElement {
            val: self.generator_val,
            field: self,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::field_element::FieldElement;

    #[test]
    fn init_test() {
        let field = galois_field!();
        assert_eq!(field.zero, 0);
        assert_eq!(field.one, 1);
        assert_eq!(field.k_modulus, 3 * u64::pow(2, 30) + 1);
        assert_eq!(field.generator_val, 5);

        assert_eq!(
            field.zero(),
            FieldElement {
                val: 0,
                field: &field
            }
        );
        assert_eq!(
            field.one(),
            FieldElement {
                val: 1,
                field: &field
            }
        );
        assert_eq!(
            field.generator(),
            FieldElement {
                val: field.generator_val,
                field: &field
            }
        );

        let test_data = [
            (-317544001, 2903681472),
            (89395626, 89395626),
            (-83021920, 3138203553),
            (-1809222732, 1412002741),
            (2884610485, 2884610485),
            (-253982795, 2967242678),
            (833431249, 833431249),
            (-1434475787, 1786749686),
            (-2604252592, 616972881),
            (1267289668, 1267289668),
            (-733379484, 2487845989),
            (-2256881083, 964344390),
            (-2368465592, 852759881),
            (249510927, 249510927),
            (-3165661992, 55563481),
            (1603141288, 1603141288),
            (-2322685483, 898539990),
            (-846661769, 2374563704),
            (-788675498, 2432549975),
            (-1101700664, 2119524809),
        ];

        for (init_val, el_val) in test_data {
            assert_eq!(field.new_element(init_val).val, el_val);
        }
    }
}
