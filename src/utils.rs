use alloc::vec::Vec;

#[inline(always)]
pub(crate) fn vec_to_vecvec<T: Clone>(
    vec: Vec<T>,
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

#[inline(always)]
pub(crate) fn vec_into<S, T>(xs: Vec<S>) -> Vec<T>
where
    S: Into<T>,
{
    xs.into_iter().map(S::into).collect()
}
