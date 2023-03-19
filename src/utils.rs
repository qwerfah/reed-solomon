use crate::{field::GaloisField, field_element::FieldElement};

use itertools::{EitherOrBoth::*, Itertools};

pub fn remove_trailing_elements<'a>(
    coeffs: &[FieldElement<'a>],
    trailing: FieldElement<'a>,
) -> Vec<FieldElement<'a>> {
    let mut filtered_coeffs = coeffs
        .iter()
        .rev()
        .skip_while(|&&el| el == trailing)
        .copied()
        .collect::<Vec<FieldElement<'a>>>();
    filtered_coeffs.reverse();
    filtered_coeffs
}

pub fn zip_longest_with_op<'a>(
    lhs: &[FieldElement<'a>],
    rhs: &[FieldElement<'a>],
    op: fn(FieldElement<'a>, FieldElement<'a>) -> FieldElement<'a>,
    fill_value: FieldElement<'a>,
) -> Vec<FieldElement<'a>> {
    lhs.iter()
        .zip_longest(rhs.iter())
        .map(|pair| match pair {
            Both(&l, &r) => op(l, r),
            Left(&l) => op(l, fill_value),
            Right(&r) => op(fill_value, r),
        })
        .collect()
}

pub fn nums_to_elements(nums: Vec<i128>, field: &GaloisField) -> Vec<FieldElement> {
    nums.into_iter().map(|num| field.new_element(num)).collect()
}
