use crate::field::GaloisField;
use std::{cmp, ops};

#[derive(Debug, Copy, Clone)]
pub struct FieldElement<'a> {
    pub val: u64,
    pub field: &'a GaloisField,
}

impl<'a> FieldElement<'a> {
    pub fn inverse(&self) -> FieldElement<'a> {
        let zero = self.field.zero as i128;
        let one = self.field.one as i128;
        let (mut t, mut new_t) = (zero, self.field.one as i128);
        let (mut r, mut new_r) = (self.field.k_modulus as i128, self.val as i128);

        while new_r != zero {
            let quotient = r / new_r;
            (t, new_t) = (new_t, t - quotient * new_t);
            (r, new_r) = (new_r, r - quotient * new_r);
        }

        assert!(r == one);
        self.field.new_element(t)
    }

    pub fn pow(self, mut n: u32) -> FieldElement<'a> {
        let mut cur_pow = self;
        let mut res = cur_pow.field.one();

        while n > 0 {
            if n % 2 != 0 {
                res *= cur_pow;
            }
            n /= 2;
            cur_pow *= cur_pow;
        }

        res
    }
}

impl<'a> ops::Add<FieldElement<'a>> for FieldElement<'a> {
    type Output = FieldElement<'a>;

    fn add(self, rhs: FieldElement<'a>) -> Self::Output {
        if std::ptr::eq(self.field, rhs.field) {
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
        self.field
            .new_element(self.field.zero as i128 - self.val as i128)
    }
}

impl<'a> ops::Sub<FieldElement<'a>> for FieldElement<'a> {
    type Output = FieldElement<'a>;

    fn sub(self, rhs: FieldElement<'a>) -> Self::Output {
        if std::ptr::eq(self.field, rhs.field) {
            self.field.new_element(self.val as i128 - rhs.val as i128)
        } else {
            panic!("Elements can't be substracted cause they lay in defferent fields");
        }
    }
}

impl<'a> ops::Mul<FieldElement<'a>> for FieldElement<'a> {
    type Output = FieldElement<'a>;

    fn mul(self, rhs: FieldElement<'a>) -> Self::Output {
        if std::ptr::eq(self.field, rhs.field) {
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
        if std::ptr::eq(self.field, rhs.field) {
            self.val = (self.val * rhs.val) % self.field.k_modulus;
        } else {
            panic!("Elements can't be multiplied cause they lay in defferent fields");
        }
    }
}

impl<'a> ops::AddAssign<FieldElement<'a>> for FieldElement<'a> {
    fn add_assign(&mut self, rhs: FieldElement<'a>) {
        if std::ptr::eq(self.field, rhs.field) {
            self.val = (self.val + rhs.val) % self.field.k_modulus;
        } else {
            panic!("Elements can't be multiplied cause they lay in defferent fields");
        }
    }
}

impl<'a> ops::SubAssign<FieldElement<'a>> for FieldElement<'a> {
    fn sub_assign(&mut self, rhs: FieldElement<'a>) {
        if std::ptr::eq(self.field, rhs.field) {
            self.val = (self.val as i128 - rhs.val as i128).rem_euclid(self.field.k_modulus as i128)
                as u64;
        } else {
            panic!("Elements can't be multiplied cause they lay in defferent fields");
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<'a> ops::Div<FieldElement<'a>> for FieldElement<'a> {
    type Output = FieldElement<'a>;

    fn div(self, rhs: FieldElement<'a>) -> Self::Output {
        if std::ptr::eq(self.field, rhs.field) {
            self * rhs.inverse()
        } else {
            panic!("Elements can't be divided cause they lay in defferent fields");
        }
    }
}

impl<'a> cmp::PartialEq<FieldElement<'a>> for FieldElement<'a> {
    fn eq(&self, other: &FieldElement<'a>) -> bool {
        std::ptr::eq(self.field, other.field) && self.val == other.val
    }
}

#[cfg(test)]
mod tests {
    use super::FieldElement;
    use crate::field::GaloisField;

    const FIELD: GaloisField = crate::galois_field!();

    #[test]
    fn inverse_test() {
        let test_data = [
            (10, 966367642),
            (60, 1234803098),
            (110, 1844883680),
            (160, 1268357530),
            (210, 2653676223),
            (260, 2019460585),
            (310, 238994148),
            (360, 2353284165),
            (410, 809234692),
            (460, 1561594088),
            (510, 2987528723),
            (560, 2203088136),
            (610, 807946717),
            (660, 2454964262),
            (710, 3189466912),
            (760, 1792866283),
            (810, 3193387722),
            (860, 1734217900),
            (910, 116813671),
            (960, 1285134746),
        ];

        for (el_val, inverse_val) in test_data {
            let element = FIELD.new_element(el_val);
            let inverse_element = FIELD.new_element(inverse_val);

            assert_eq!(element.inverse(), inverse_element);
        }
    }

    #[test]
    fn pow_test() {
        let test_data = [
            (-2437383495, 81, 3187458499),
            (-2661509949, 88, 999582370),
            (2373853268, 79, 843400527),
            (-1712893553, 80, 2630252157),
            (-1056252449, 29, 852554515),
            (2316116930, 22, 1942594479),
            (1039791324, 20, 721898573),
            (2110495054, 100, 1382847113),
            (-176065273, 77, 1658860552),
            (1458457238, 49, 716007232),
            (-697405797, 21, 83985176),
            (2563556474, 14, 2031213162),
            (-2084353984, 91, 2599022742),
            (32691010, 94, 2102602746),
            (-982274076, 25, 2868041609),
            (-1504968126, 80, 1216473615),
            (-2959700505, 44, 505173613),
            (-387879001, 60, 923994655),
            (903348446, 81, 1007376090),
            (-1652511619, 42, 2392825049),
        ];

        for (el_val, pow, pow_val) in test_data {
            assert_eq!(FIELD.new_element(el_val).pow(pow).val, pow_val);
        }
    }

    #[test]
    fn sum_test() {
        let test_data = [
            (-2674122163, -2977635245, 790693538),
            (2576410027, 2304397362, 1659581916),
            (951208061, 1011963981, 1963172042),
            (2226545896, -620302777, 1606243119),
            (-2557274478, -3157155296, 728021172),
            (-851027498, 158235079, 2528433054),
            (639039827, 2867198104, 285012458),
            (3130882677, -1285797189, 1845085488),
            (1951540482, 2991978173, 1722293182),
            (-1918350070, -1984437146, 2539663730),
            (-340934393, -319180320, 2561110760),
            (-1501705013, -625605485, 1093914975),
            (487475469, -1652566431, 2056134511),
            (1571979455, -313317621, 1258661834),
            (2800761065, 1742772228, 1322307820),
            (365719652, -1689462048, 1897483077),
            (-2364878684, 152486048, 1008832837),
            (-2973372868, -2269352785, 1199725293),
            (-732571295, 2695959969, 1963388674),
            (1029418390, 2673697386, 481890303),
        ];

        for (lhs_val, rhs_val, sum_val) in test_data {
            let lhs = FIELD.new_element(lhs_val);
            let rhs = FIELD.new_element(rhs_val);
            assert_eq!((lhs + rhs).val, sum_val);
        }
    }

    #[test]
    fn sub_test() {
        let test_data = [
            (-1086477693, -2867903489, 1781425796),
            (1462968467, 2104794817, 2579399123),
            (-981229042, 711386736, 1528609695),
            (-3122116852, 1019038246, 2301295848),
            (238649413, -2006128286, 2244777699),
            (224553338, -1631995338, 1856548676),
            (-2902747518, 2222784665, 1316918763),
            (-1970920006, 643853230, 606452237),
            (3182612939, 1092239114, 2090373825),
            (-2629398192, -2454056762, 3045884043),
            (281323751, 1791809146, 1710740078),
            (-606349471, -1943283012, 1336933541),
            (683948720, -882251417, 1566200137),
            (2113760060, -3006792663, 1899327250),
            (-1892589503, -2215774910, 323185407),
            (-2759014685, 1531742391, 2151693870),
            (2138384011, 3187835787, 2171773697),
            (-1414458833, 1861660164, 3166331949),
            (-237701752, 1786469979, 1197053742),
            (-186446840, 545421759, 2489356874),
        ];

        for (lhs_val, rhs_val, sub_val) in test_data {
            let lhs = FIELD.new_element(lhs_val);
            let rhs = FIELD.new_element(rhs_val);
            assert_eq!((lhs - rhs).val, sub_val);
        }
    }

    #[test]
    fn neg_test() {
        let test_data = [
            (513023168, 2708202305),
            (-1567894667, 1567894667),
            (-351452192, 351452192),
            (490763083, 2730462390),
            (2477566028, 743659445),
            (-788183888, 788183888),
            (1846613209, 1374612264),
            (-1718246028, 1718246028),
            (-2286428665, 2286428665),
            (-246726409, 246726409),
            (591661579, 2629563894),
            (2659897702, 561327771),
            (-2959779295, 2959779295),
            (2679448172, 541777301),
            (-1464856126, 1464856126),
            (2011683406, 1209542067),
            (288395027, 2932830446),
            (554377246, 2666848227),
            (-1610816554, 1610816554),
            (431477634, 2789747839),
        ];

        for (el_val, neg_val) in test_data {
            assert_eq!((-FIELD.new_element(el_val)).val, neg_val);
        }
    }

    #[test]
    fn mul_test() {
        let test_data = [
            (-3217314898, 1575806093, 503069339),
            (2848844329, 2969028563, 1187663836),
            (3062075925, 2178575679, 788913012),
            (-700399855, -990740140, 1119959245),
            (-2363483117, 264996178, 2055646916),
            (-1346656139, 2643135521, 1836704055),
            (-1393117006, -1577764739, 2022522856),
            (2794778151, -1991645734, 2197392325),
            (-716994424, 1497407377, 3135933462),
            (-1267213260, -2282175351, 1637441737),
            (-2875408142, -1537917493, 2026713801),
            (-2197168149, 466380145, 2716701938),
            (-580535173, 2820609132, 47625444),
            (-322134490, 1368935197, 3124491378),
            (2445947494, 2445970384, 426672109),
            (-2631271225, 1089981044, 1770934086),
            (-1519254363, 1676274862, 2282999215),
            (254581814, 787315415, 2402383997),
            (2533391896, -1240172834, 1351560160),
            (1097669165, -794986387, 1836571338),
        ];

        for (lhs_val, rhs_val, mul_val) in test_data {
            let lhs = FIELD.new_element(lhs_val);
            let rhs = FIELD.new_element(rhs_val);
            assert_eq!((lhs * rhs).val, mul_val);
        }
    }

    #[test]
    fn div_test() {
        let test_data = [
            (-408821234, -2199589188, 811607262),
            (2160990732, 2051454364, 727573819),
            (-2047046118, -2526017147, 944282317),
            (459438645, -2727940221, 876871549),
            (-2535264262, -1022515841, 89145159),
            (-2583698943, 1561204688, 1283822742),
            (766084880, 1324295007, 2980830939),
            (1766791713, -1409952932, 309991435),
            (-121584279, 2439279458, 2619494963),
            (2842321123, 1547202657, 851510893),
            (1432991782, -112282547, 474778458),
            (740034210, -1605762111, 458854688),
            (1940358425, -1199047992, 348921429),
            (2648259158, -2602894318, 2123343152),
            (393609078, -1782376871, 519108401),
            (-1435717664, 2101411704, 903158850),
            (-1149589941, 2835367238, 2875945717),
            (-1092076334, 791436329, 1001337653),
            (2556250327, -1448960692, 3109683912),
            (94489230, -1256612031, 3104814716),
        ];

        for (lhs_val, rhs_val, div_val) in test_data {
            let lhs = FIELD.new_element(lhs_val);
            let rhs = FIELD.new_element(rhs_val);
            assert_eq!((lhs / rhs).val, div_val);
        }
    }

    #[test]
    fn expr_test() {
        let element = FIELD.new_element(2).pow(30) * FIELD.new_element(3) + FIELD.new_element(1);
        assert_eq!(
            element,
            FieldElement {
                val: 0,
                field: &FIELD
            }
        );
    }
}
