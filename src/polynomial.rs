use crate::field::GaloisField;
use crate::field_element::FieldElement;
use crate::utils;

use std::cmp;
use std::ops;

/// Polynomial above some finite field `field`.
/// # Arguments
/// * `coeffs` - the coefficients of the polynomial, listed
/// in ascending order of terms powers (they all should lay in the same field `field`)
/// * `field` - some finite field that provides polynomial coefficients
/// * `var` - polynomial variable designation
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
            coeffs: utils::remove_trailing_elements(&coeffs, zero),
            field,
            var: var.to_string(),
        }
    }

    pub fn from(coeffs: Vec<FieldElement<'a>>, other: &'a Polynomial) -> Polynomial<'a> {
        Polynomial::new(coeffs, other.field, &other.var)
    }

    pub fn empty(other: &'a Polynomial) -> Polynomial<'a> {
        Polynomial::new(vec![], other.field, &other.var)
    }

    pub fn get_field_ptr(coeffs: &[FieldElement]) -> *const GaloisField {
        match coeffs.first() {
            Some(element) => element.field as *const _,
            _ => std::ptr::null() as *const _,
        }
    }

    pub fn deg(&self) -> i64 {
        self.coeffs.len() as i64 - 1
    }

    pub fn non_empty(&self) -> bool {
        !self.coeffs.is_empty()
    }

    /// Function composition operation on two polynomials.
    pub fn compose(&self, rhs: Polynomial<'a>) -> Polynomial {
        let mut res = Polynomial::empty(self);

        for &coef in self.coeffs.iter().rev() {
            res = (res * rhs.clone()) + Polynomial::new(vec![coef], self.field, &self.var);
        }

        res
    }

    pub fn qdiv(&self, rhs: &Polynomial<'a>) -> (Polynomial, Polynomial) {
        Polynomial::check_bin_op_args(self, rhs);
        self.qdiv_(rhs)
    }

    pub fn monomial(deg: usize, coef: FieldElement<'a>, field: &'a GaloisField) -> Polynomial<'a> {
        let mut coeffs = vec![field.zero(); deg];
        coeffs.push(coef);
        Polynomial::new(coeffs, field, "x")
    }

    pub fn interpolate(
        x: &'a Vec<FieldElement<'a>>,
        y: &'a Vec<FieldElement<'a>>,
    ) -> Polynomial<'a> {
        if x.is_empty() || y.is_empty() || x.len() != y.len() {
            panic!("Impossible to establish a one-to-one correspondence between the definition and value areas.");
        }

        let expected_field_ptr = x.first().unwrap().field as *const _;

        for elem in x.iter().chain(y.iter()) {
            if elem.field as *const _ != expected_field_ptr {
                panic!("Not all elements are lay in the same field.");
            }
        }

        let polynomials = Polynomial::calculate_lagrange_polynomials(x);
        Polynomial::interpolate_poly_lagrange(y, polynomials)
    }

    #[allow(unused)]
    fn calculate_lagrange_polynomials(x: &'a [FieldElement<'a>]) -> Vec<Polynomial<'a>> {
        unimplemented!();
    }

    #[allow(unused)]
    fn interpolate_poly_lagrange(
        y: &'a [FieldElement<'a>],
        polynomials: Vec<Polynomial<'a>>,
    ) -> Polynomial<'a> {
        unimplemented!();
    }

    /// Calculates quotient and remainder polynomials such that
    /// f = q * g + r, where deg(r) < deg(g).
    fn qdiv_(&self, rhs: &Polynomial<'a>) -> (Polynomial, Polynomial) {
        let rhs_coeffs = utils::remove_trailing_elements(&rhs.coeffs, rhs.field.zero());
        assert!(!rhs_coeffs.is_empty());

        let lhs_coeffs = utils::remove_trailing_elements(&self.coeffs, self.field.zero());
        if lhs_coeffs.is_empty() {
            (Polynomial::empty(self), Polynomial::empty(self))
        } else {
            let mut rem = lhs_coeffs;
            let mut deg_dif = rem.len() as i128 - rhs_coeffs.len() as i128;
            let quotient_len = if deg_dif >= -1 {
                (deg_dif + 1) as usize
            } else {
                0
            };
            let mut quotient = vec![self.field.zero(); quotient_len];
            let g_msc_inv = rhs_coeffs.last().unwrap().inverse();

            while deg_dif >= 0 {
                let tmp = *rem.last().unwrap() * g_msc_inv;
                let mut last_non_zero = deg_dif - 1;

                quotient[deg_dif as usize] += tmp;

                let offset = deg_dif as usize;
                for (i, coef) in rhs_coeffs.iter().enumerate() {
                    rem[i + offset] -= tmp * *coef;
                    if rem[i + offset] != self.field.zero() {
                        last_non_zero = (i + offset) as i128;
                    }
                }

                rem.truncate((last_non_zero + 1) as usize);
                deg_dif = rem.len() as i128 - rhs_coeffs.len() as i128;
            }

            (
                Polynomial::from(
                    utils::remove_trailing_elements(&quotient, self.field.zero()),
                    self,
                ),
                Polynomial::from(rem, self),
            )
        }
    }

    fn check_bin_op_args(lhs: &Polynomial, rhs: &Polynomial) {
        if lhs.field as *const _ != rhs.field as *const _ {
            panic!("Polynomials are biult over different fields!");
        }

        if lhs.var != rhs.var {
            panic!("Polynomials have different variable names!");
        }
    }

    fn bin_op(
        lhs: &Polynomial<'a>,
        rhs: &Polynomial<'a>,
        op: fn(FieldElement<'a>, FieldElement<'a>) -> FieldElement<'a>,
    ) -> Polynomial<'a> {
        Polynomial::check_bin_op_args(lhs, rhs);

        Polynomial {
            coeffs: utils::zip_longest_with_op(&lhs.coeffs, &rhs.coeffs, op, lhs.field.zero()),
            field: lhs.field,
            var: lhs.var.to_string(),
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
        Polynomial::bin_op(&self, &rhs, |a, b| a + b)
    }
}

impl<'a> ops::Sub<Polynomial<'a>> for Polynomial<'a> {
    type Output = Polynomial<'a>;

    fn sub(self, rhs: Polynomial<'a>) -> Self::Output {
        Polynomial::bin_op(&self, &rhs, |a, b| a - b)
    }
}

impl<'a> ops::Neg for Polynomial<'a> {
    type Output = Polynomial<'a>;

    fn neg(self) -> Self::Output {
        Polynomial::new(vec![], self.field, &self.var) - self
    }
}

impl<'a> ops::Mul<Polynomial<'a>> for Polynomial<'a> {
    type Output = Polynomial<'a>;

    fn mul(self, rhs: Polynomial<'a>) -> Self::Output {
        Polynomial::check_bin_op_args(&self, &rhs);

        let lhs_raw_coeffs = self
            .coeffs
            .iter()
            .map(|elem| elem.val)
            .collect::<Vec<u64>>();
        let rhs_raw_coeffs = rhs.coeffs.iter().map(|elem| elem.val).collect::<Vec<u64>>();
        let res_len = self.deg() + rhs.deg() + 1;
        let mut res_raw_coeffs = vec![0; cmp::max(res_len, 0) as usize];

        for (i, lhs_val) in lhs_raw_coeffs.into_iter().enumerate() {
            for (j, rhs_val) in rhs_raw_coeffs.iter().enumerate() {
                res_raw_coeffs[i + j] += lhs_val * rhs_val;
            }
        }

        Polynomial::new(
            res_raw_coeffs
                .into_iter()
                .map(|val| self.field.new_element(val as i128))
                .collect(),
            self.field,
            &self.var,
        )
    }
}

impl<'a> ops::Div<&Polynomial<'a>> for &'a Polynomial<'a> {
    type Output = Polynomial<'a>;

    fn div(self, rhs: &Polynomial<'a>) -> Self::Output {
        Polynomial::check_bin_op_args(self, rhs);
        let (quot, rem) = self.qdiv_(rhs);
        if rem.coeffs.is_empty() {
            quot
        } else {
            panic!("Polynomials are not divisible.");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Polynomial;
    use crate::field::GaloisField;
    use crate::utils;

    use itertools::izip;

    const FIELD: GaloisField = crate::galois_field!();

    #[test]
    fn init_test() {
        assert_eq!(
            Polynomial::x(&FIELD),
            Polynomial {
                coeffs: vec![FIELD.zero(), FIELD.one()],
                field: &FIELD,
                var: "x".to_string()
            }
        );
    }

    #[test]
    fn add_test() {
        let empty_poly = Polynomial::new(vec![], &FIELD, "x");
        let some_poly = Polynomial::new(
            vec![
                FIELD.new_element(1),
                FIELD.new_element(2),
                FIELD.new_element(3),
            ],
            &FIELD,
            "x",
        );
        let other_poly = Polynomial::new(
            vec![
                FIELD.new_element(2),
                FIELD.new_element(3),
                FIELD.new_element(4),
                FIELD.new_element(5),
                FIELD.new_element(6),
            ],
            &FIELD,
            "x",
        );

        let sum1 = Polynomial::new(
            vec![
                FIELD.new_element(1),
                FIELD.new_element(2),
                FIELD.new_element(3),
            ],
            &FIELD,
            "x",
        );

        let sum2 = Polynomial::new(
            vec![
                FIELD.new_element(3),
                FIELD.new_element(5),
                FIELD.new_element(7),
                FIELD.new_element(5),
                FIELD.new_element(6),
            ],
            &FIELD,
            "x",
        );

        assert_eq!(empty_poly.clone() + empty_poly.clone(), empty_poly);
        assert_eq!(empty_poly.clone() + some_poly.clone(), sum1);
        assert_eq!(some_poly.clone() + empty_poly, sum1);
        assert_eq!(some_poly.clone() + other_poly.clone(), sum2);
        assert_eq!(other_poly + some_poly, sum2);
    }

    #[test]
    fn sub_test() {
        let lhs_data = [
            vec![129, 146, 100, 132, 182, 174, 153, 185, 147],
            vec![152, 100, 128, 138],
            vec![110, 142, 143, 109, 184, 147],
            vec![114, 137, 199, 170, 161],
            vec![181, 136, 160, 179, 188, 144, 141, 137, 121],
            vec![118, 186, 146, 173, 139, 125],
            vec![182, 115, 147, 183, 116, 161, 193, 163, 153],
            vec![101, 158, 185, 164, 170, 111, 118, 189, 181],
            vec![186, 125, 158, 116],
            vec![146, 116, 149, 143, 200, 149],
        ];
        let rhs_data = [
            vec![144, 134, 169, 160, 153, 160],
            vec![119, 183, 157, 118, 134],
            vec![141, 150, 200, 131, 148, 144, 193, 118],
            vec![187, 126, 176, 150, 106, 176, 187, 138, 102, 102],
            vec![148, 141, 141, 180, 106, 128, 134],
            vec![165, 117, 153, 121, 150, 152, 134, 155, 181],
            vec![155, 146, 171],
            vec![101, 106, 180, 111, 114, 188, 103, 144],
            vec![154, 134, 142, 119, 139, 152, 140, 177, 166],
            vec![132, 164, 173, 171, 198, 102, 198, 196, 194],
        ];
        let res_data = [
            vec![-15, 12, -69, -28, 29, 14, 153, 185, 147],
            vec![33, -83, -29, 20, -134],
            vec![-31, -8, -57, -22, 36, 3, -193, -118],
            vec![-73, 11, 23, 20, 55, -176, -187, -138, -102, -102],
            vec![33, -5, 19, -1, 82, 16, 7, 137, 121],
            vec![-47, 69, -7, 52, -11, -27, -134, -155, -181],
            vec![27, -31, -24, 183, 116, 161, 193, 163, 153],
            vec![0, 52, 5, 53, 56, -77, 15, 45, 181],
            vec![32, -9, 16, -3, -139, -152, -140, -177, -166],
            vec![14, -48, -24, -28, 2, 47, -198, -196, -194],
        ];

        for (lhs_raw, rhs_raw, res_raw) in izip!(lhs_data, rhs_data, res_data) {
            let (lhs_poly, rhs_poly, res_poly) = prepare_data_for_bin_op(lhs_raw, rhs_raw, res_raw);
            assert_eq!(lhs_poly - rhs_poly, res_poly);
        }
    }

    #[test]
    fn neg_test() {
        let poly_data = [
            vec![180, 144, 113, 100, 172],
            vec![188, 108, 152, 151],
            vec![121, 171, 195, 160, 157],
            vec![151, 156, 182, 168, 123, 104, 148],
            vec![138, 181, 123, 106, 157, 116, 140, 127, 147],
            vec![102, 186, 130, 115],
            vec![155, 101, 116, 162],
            vec![167, 115, 115, 185],
            vec![163, 112, 184, 152, 126, 185, 198, 198, 164, 175],
            vec![178, 132, 157, 155],
        ];
        let res_data = [
            vec![-180, -144, -113, -100, -172],
            vec![-188, -108, -152, -151],
            vec![-121, -171, -195, -160, -157],
            vec![-151, -156, -182, -168, -123, -104, -148],
            vec![-138, -181, -123, -106, -157, -116, -140, -127, -147],
            vec![-102, -186, -130, -115],
            vec![-155, -101, -116, -162],
            vec![-167, -115, -115, -185],
            vec![-163, -112, -184, -152, -126, -185, -198, -198, -164, -175],
            vec![-178, -132, -157, -155],
        ];

        for (poly_raw, res_raw) in izip!(poly_data, res_data) {
            let poly_coeffs = utils::nums_to_elements(poly_raw, &FIELD);
            let res_coeffs = utils::nums_to_elements(res_raw, &FIELD);

            let poly = Polynomial::new(poly_coeffs, &FIELD, "x");
            let res = Polynomial::new(res_coeffs, &FIELD, "x");

            assert_eq!(-poly, res);
        }
    }

    #[test]
    fn mul_test() {
        let lhs_data = [
            vec![16, 10, 13],
            vec![12, 18],
            vec![19, 20, 12, 10, 15],
            vec![12, 20, 18, 17],
            vec![13, 12, 12, 19],
            vec![16, 11],
            vec![11, 18, 11, 12, 16],
            vec![14, 20, 11, 19],
            vec![15, 12, 18, 18],
            vec![14, 12, 13, 13],
        ];
        let rhs_data = [
            vec![12, 18, 20, 15],
            vec![11, 14],
            vec![19, 16, 10, 18, 10],
            vec![19, 20, 10],
            vec![19, 18],
            vec![17, 13, 20],
            vec![11, 12, 17, 18],
            vec![18, 15, 18, 14],
            vec![16, 12, 13, 12],
            vec![10, 13, 10, 12],
        ];
        let res_data = [
            vec![192, 408, 656, 674, 410, 195],
            vec![132, 366, 252],
            vec![361, 684, 738, 924, 1115, 756, 450, 370, 150],
            vec![228, 620, 862, 883, 520, 170],
            vec![247, 462, 444, 577, 342],
            vec![272, 395, 463, 220],
            vec![121, 330, 524, 768, 831, 594, 488, 288],
            vec![252, 570, 750, 1063, 763, 496, 266],
            vec![240, 372, 627, 840, 594, 450, 216],
            vec![140, 302, 426, 587, 443, 286, 156],
        ];

        for (lhs_raw, rhs_raw, res_raw) in izip!(lhs_data, rhs_data, res_data) {
            let (lhs_poly, rhs_poly, res_poly) = prepare_data_for_bin_op(lhs_raw, rhs_raw, res_raw);
            assert_eq!(lhs_poly * rhs_poly, res_poly);
        }
    }

    #[test]
    fn compose_test() {
        let lhs_data = [
            vec![17, 10, 12, 19, 19],
            vec![20, 13],
            vec![18, 13, 10, 15, 11],
            vec![18, 11],
            vec![10, 18],
            vec![18, 20, 19],
            vec![19, 17, 10, 15],
            vec![13, 10, 16, 19, 13],
            vec![19, 18],
            vec![17, 15],
        ];
        let rhs_data = [
            vec![20, 13],
            vec![11, 18],
            vec![20, 20, 18],
            vec![13, 11, 12, 16, 15],
            vec![15, 16, 14],
            vec![16, 15, 16],
            vec![12, 14, 13, 16, 12],
            vec![17, 16, 18],
            vec![13, 20, 12, 14, 17],
            vec![19, 11, 14, 10],
        ];
        let res_data = vec![
            vec![3197017, 8206770, 7901088, 3381183, 542659],
            vec![163, 234],
            vec![
                1884278, 7408260, 17591434, 26823200, 29940440, 23734800, 13773240, 5132160,
                1154736,
            ],
            vec![161, 121, 132, 176, 165],
            vec![280, 288, 252],
            vec![5202, 9420, 14323, 9120, 4864],
            vec![
                27583, 94318, 195381, 349152, 534854, 661190, 713275, 690720, 562140, 376800,
                222480, 103680, 25920,
            ],
            vec![
                1183927, 4360048, 10927966, 17250240, 20957788, 17877312, 11735928, 4852224,
                1364688,
            ],
            vec![253, 360, 216, 252, 306],
            vec![302, 165, 210, 150],
        ];

        for (lhs_raw, rhs_raw, res_raw) in izip!(lhs_data, rhs_data, res_data) {
            let (lhs_poly, rhs_poly, res_poly) = prepare_data_for_bin_op(lhs_raw, rhs_raw, res_raw);
            assert_eq!(lhs_poly.compose(rhs_poly), res_poly);
        }
    }

    #[test]
    fn qdiv_test() {
        let lhs_data = [
            vec![176, 157, 149, 103],
            vec![177, 166, 195, 122, 115, 145],
            vec![117, 110, 127],
            vec![102, 192, 116, 139],
            vec![169, 103, 195, 127, 141],
            vec![108, 142, 150],
            vec![119, 138, 195, 108, 127, 132],
            vec![113, 129, 128, 111, 160, 139],
            vec![191, 146, 126, 186, 106],
            vec![170, 182, 141, 138, 180, 131, 129],
        ];
        let rhs_data = [
            vec![123, 184, 191, 100],
            vec![129, 127, 140, 133, 194, 100, 169],
            vec![136, 124, 147],
            vec![151, 173, 158, 194, 178, 132],
            vec![164, 107, 122, 168, 103, 162],
            vec![145, 193, 178, 101, 117, 177],
            vec![185, 139, 116],
            vec![135, 191, 107, 133, 129, 120, 106],
            vec![158, 138, 154, 132, 197],
            vec![139, 144, 167, 171, 116, 155],
        ];
        let res_data = [
            (vec![-354334801], vec![-1513975923, 773094081, 32212207]),
            (vec![], vec![177, 166, 195, 122, 115, 145]),
            (vec![-810784642], vec![745045347, 679306055]),
            (vec![], vec![102, 192, 116, 139]),
            (vec![], vec![169, 103, 195, 127, 141]),
            (vec![], vec![108, 142, 150]),
            (
                vec![153572227, -394679037, -95755811, -777537182],
                vec![580167381, 129474862],
            ),
            (vec![], vec![113, 129, 128, 111, 160, 139]),
            (
                vec![310676569],
                vec![-768515616, -997435227, 474190595, 866624227],
            ),
            (
                vec![426770476, 457206197],
                vec![-1339037480, 621183702, 1404333620, -1154814399, 1161384197],
            ),
        ];

        for (lhs_raw, rhs_raw, (quot_raw, rem_raw)) in izip!(lhs_data, rhs_data, res_data) {
            let (lhs_poly, rhs_poly, quot_poly) =
                prepare_data_for_bin_op(lhs_raw, rhs_raw, quot_raw);
            let rem_coeffs = utils::nums_to_elements(rem_raw, &FIELD);
            let rem_poly = Polynomial::new(rem_coeffs, &FIELD, "x");

            assert_eq!(lhs_poly.qdiv(&rhs_poly), (quot_poly, rem_poly));
        }
    }

    #[test]
    fn monomial_test() {
        assert_eq!(
            Polynomial::monomial(0, FIELD.new_element(12), &FIELD).coeffs,
            vec![FIELD.new_element(12)]
        );
        assert_eq!(
            Polynomial::monomial(3, FIELD.new_element(12), &FIELD).coeffs,
            vec![
                FIELD.new_element(0),
                FIELD.new_element(0),
                FIELD.new_element(0),
                FIELD.new_element(12)
            ]
        );
    }

    fn prepare_data_for_bin_op<'a>(
        lhs_raw: Vec<i128>,
        rhs_raw: Vec<i128>,
        res_raw: Vec<i128>,
    ) -> (Polynomial<'a>, Polynomial<'a>, Polynomial<'a>) {
        let lhs_coeffs = utils::nums_to_elements(lhs_raw, &FIELD);
        let rhs_coeffs = utils::nums_to_elements(rhs_raw, &FIELD);
        let res_coeffs = utils::nums_to_elements(res_raw, &FIELD);
        let lhs_poly = Polynomial::new(lhs_coeffs, &FIELD, "x");
        let rhs_poly = Polynomial::new(rhs_coeffs, &FIELD, "x");
        let res_poly = Polynomial::new(res_coeffs, &FIELD, "x");

        (lhs_poly, rhs_poly, res_poly)
    }
}
