use crate::field_element::FieldElement;

use itertools::{EitherOrBoth::*, Itertools};

pub fn remove_trailing_elements<'a>(
    coeffs: Vec<FieldElement<'a>>,
    trailing: FieldElement<'a>,
) -> Vec<FieldElement<'a>> {
    let mut filtered_coeffs = coeffs
        .into_iter()
        .rev()
        .skip_while(|el| *el == trailing)
        .collect::<Vec<FieldElement<'a>>>();
    filtered_coeffs.reverse();
    filtered_coeffs
}

pub fn zip_longest_with_op<'a>(
    lhs: &Vec<FieldElement<'a>>,
    rhs: &Vec<FieldElement<'a>>,
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
