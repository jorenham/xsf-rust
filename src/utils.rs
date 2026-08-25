#[inline]
pub(crate) fn out2<T: Copy + From<f64>>(f: impl FnOnce(*mut T, *mut T)) -> (T, T) {
    let (mut a, mut b) = (T::from(f64::NAN), T::from(f64::NAN));
    f(&raw mut a, &raw mut b);
    (a, b)
}

#[inline]
pub(crate) fn out4<T: Copy + From<f64>>(
    f: impl FnOnce(*mut T, *mut T, *mut T, *mut T),
) -> (T, T, T, T) {
    let nan = T::from(f64::NAN);
    let (mut a, mut b, mut c, mut d) = (nan, nan, nan, nan);
    f(&raw mut a, &raw mut b, &raw mut c, &raw mut d);
    (a, b, c, d)
}

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

#[inline]
pub(crate) fn vec_into_vecvec<S: Into<T>, T: Clone>(
    vec: Vec<S>,
    rows: usize,
    cols: usize,
    transpose: bool,
) -> Vec<Vec<T>> {
    vec_to_vecvec(&vec_into(vec), rows, cols, transpose)
}
