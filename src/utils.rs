/// ln(pi)
pub(crate) const LN_PI: f64 = 1.144_729_885_849_400_2;

#[inline]
pub(crate) fn vec_into<S: Into<T>, T>(xs: Vec<S>) -> Vec<T> {
    xs.into_iter().map(S::into).collect()
}

#[inline]
pub(crate) fn vec_to_vecvec<T: Clone>(
    vec: &[T],
    rows: usize,
    cols: usize,
    transpose: bool,
) -> Vec<Vec<T>> {
    if transpose {
        // Transpose: convert from (rows, cols) to (cols, rows)
        (0..cols)
            .map(|j| (0..rows).map(|i| vec[i * cols + j].clone()).collect())
            .collect()
    } else {
        // Normal: convert from flat vec to (rows, cols)
        (0..rows)
            .map(|i| vec[i * cols..(i + 1) * cols].to_vec())
            .collect()
    }
}

/// Keep the first `rows` rows and `cols` columns
#[inline]
pub(crate) fn vecvec_truncate<T>(mut vv: Vec<Vec<T>>, rows: usize, cols: usize) -> Vec<Vec<T>> {
    vv.truncate(rows);
    for v in &mut vv {
        v.truncate(cols);
    }
    vv
}

#[inline]
pub(crate) fn vec_into_vecvec<S: Into<T>, T: Clone>(
    vec: Vec<S>,
    rows: usize,
    cols: usize,
    transpose: bool,
) -> Vec<Vec<T>> {
    vec_to_vecvec(&vec_into(vec), rows, cols, transpose)
}
