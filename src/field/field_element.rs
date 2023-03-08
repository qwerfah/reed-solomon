use crate::field::field::GaloisField;
use std::ops;

#[derive(Debug, Clone)]
pub struct FieldElement<'a> {
    pub val: u32,
    pub field: &'a GaloisField,
}

impl<'a> FieldElement<'a> {
    pub fn inverse(self) -> FieldElement<'a> {
        let (mut t, mut new_t) = (self.field.zero, self.field.one);
        let (mut r, mut new_r) = (self.field.k_modulus, self.val);

        while new_r != self.field.zero {
            let quotient = r / new_r;
            (t, new_t) = (new_t, (t - (quotient * new_t)));
            (r, new_r) = (new_r, r - quotient * new_r);
        }
        assert!(r == self.field.one);

        FieldElement {
            val: t,
            field: self.field,
        }
    }

    pub fn pow(self, mut n: u32) -> FieldElement<'a> {
        let mut cur_pow = self;
        let mut res = cur_pow.field.one();

        while n > 0 {
            if n % 2 != 0 {
                res = res * cur_pow;
            }
            n /= 2;
            cur_pow = cur_pow.clone() * cur_pow.clone();
        }

        res
    }
}

impl<'a> ops::Add<FieldElement<'a>> for FieldElement<'a> {
    type Output = FieldElement<'a>;

    fn add(self, rhs: FieldElement<'a>) -> Self::Output {
        if self.field as *const _ == rhs.field as *const _ {
            FieldElement {
                val: (self.val + rhs.val) % self.field.k_modulus,
                field: self.field,
            }
        } else {
            panic!("Elements can't be summed cause they lay in defferent fields");
        }
    }
}

impl<'a> ops::Neg for FieldElement<'a> {
    type Output = FieldElement<'a>;

    fn neg(self) -> FieldElement<'a> {
        FieldElement {
            val: (self.field.zero - self.val) % self.field.k_modulus,
            field: self.field,
        }
    }
}

impl<'a> ops::Sub<FieldElement<'a>> for FieldElement<'a> {
    type Output = FieldElement<'a>;

    fn sub(self, rhs: FieldElement<'a>) -> Self::Output {
        if self.field as *const _ == rhs.field as *const _ {
            FieldElement {
                val: (self.val - rhs.val) % self.field.k_modulus,
                field: self.field,
            }
        } else {
            panic!("Elements can't be substracted cause they lay in defferent fields");
        }
    }
}

impl<'a> ops::Mul<FieldElement<'a>> for FieldElement<'a> {
    type Output = FieldElement<'a>;

    fn mul(self, rhs: FieldElement<'a>) -> Self::Output {
        if self.field as *const _ == rhs.field as *const _ {
            FieldElement {
                val: (self.val * rhs.val) % self.field.k_modulus,
                field: self.field,
            }
        } else {
            panic!("Elements can't be multiplied cause they lay in defferent fields");
        }
    }
}

impl<'a> ops::MulAssign<FieldElement<'a>> for FieldElement<'a> {
    fn mul_assign(&mut self, rhs: FieldElement<'a>) {
        if self.field as *const _ == rhs.field as *const _ {
            self.val = (self.val * rhs.val) % self.field.k_modulus;
        } else {
            panic!("Elements can't be multiplied cause they lay in defferent fields");
        }
    }
}

impl<'a> ops::Div<FieldElement<'a>> for FieldElement<'a> {
    type Output = FieldElement<'a>;

    fn div(self, rhs: FieldElement<'a>) -> Self::Output {
        if self.field as *const _ == rhs.field as *const _ {
            self * rhs.inverse()
        } else {
            panic!("Elements can't be divided cause they lay in defferent fields");
        }
    }
}
