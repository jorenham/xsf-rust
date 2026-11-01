use num_traits::{PrimInt, Unsigned};

pub trait Stirling2Output: crate::sealed::Sealed {
    fn stirling2(n: u32, k: u32) -> Self;
}

macro_rules! impl_stirling2_exact {
    ($($T:ty),*) => {$(
        impl Stirling2Output for $T {
            #[inline]
            fn stirling2(n: u32, k: u32) -> Self {
                stirling2_exact(n, k)
            }
        }
    )*};
}

impl_stirling2_exact!(u8, u16, u32, u64, u128, usize);

impl Stirling2Output for f64 {
    #[inline]
    fn stirling2(n: u32, k: u32) -> Self {
        crate::xsf::stirling2_float(n.into(), k.into())
    }
}

/// Stirling number of the second kind *S(n,k)*
///
/// Stirling numbers of the second kind count the number of ways to partition a set with *n*
/// elements into *k* non-empty subsets.
///
/// For an unsigned integer return type `T`, this is a pure rust translation of
/// [`scipy.special.stirling2(n, k, exact=True)`][stirling2].
///
/// For `T = f64`, this corresponds to [`scipy.special.stirling2(n, k, exact=False)`][stirling2],
/// which for *n <= 50* uses dynamic programming, and for *n > 50* the second-order asymptotic
/// approximation of Temme \[1\]. The latter has a relative error of up to about `2e-5` (at
/// *n = 51*).
///
/// [stirling2]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stirling2.html
///
/// \[1\] N. M. Temme, "Asymptotic estimates of Stirling numbers", Studies in Applied Mathematics,
///     89(3):233-243, 1993.
///
/// # Examples
///
/// ```
/// use xsf::stirling2;
///
/// assert_eq!(stirling2::<u64>(10, 3), 9_330);
/// assert_eq!(stirling2::<u64>(10, 0), 0);
/// assert_eq!(stirling2::<u64>(9, 3), 3_025);
///
/// // won't fit in u64
/// assert_eq!(stirling2::<u128>(42, 4), 805_932_309_912_663_709_372_025);
///
/// // won't fit in u128
/// assert_eq!(stirling2::<f64>(50, 10), 2.615_471_651_586_288e43);
/// ```
///
/// # Types
/// - `T` (required): unsigned integer or `f64` return type
///
/// # See also
/// - [`comb`](crate::comb): *k*-combinations of *n* things, <sub>*n*</sub>C<sub>*k*</sub>
/// - [`perm`](crate::perm): *k*-permutations of *n* things, <sub>*n*</sub>P<sub>*k*</sub>
///
/// # Panics
/// - For an integer `T`, this function will panic if the result does not fit in type `T`.
#[must_use]
#[inline]
pub fn stirling2<T: Stirling2Output>(n: u32, k: u32) -> T {
    T::stirling2(n, k)
}

#[inline]
fn stirling2_exact<T: PrimInt + Unsigned>(n: u32, k: u32) -> T {
    if n == k || (n > 0 && k == 1) {
        return T::one();
    } else if n < k || k == 0 {
        return T::zero();
    } else if k == 2 {
        // S(n, 2) = 2^(n-1) - 1
        return 129_u32
            .checked_sub(n)
            .and_then(|shift| T::from(u128::MAX >> shift))
            .unwrap();
    } else if n == k + 1 {
        // S(n, n-1) = C(n, 2) = n(n-1)/2
        return T::from(u128::from(n) * (u128::from(n) - 1) / 2).unwrap();
    }

    let k = k as usize;
    let mut n_row = vec![T::one(); k];
    for i in 2..n as usize {
        for j in (1..i.min(k)).rev() {
            n_row[j] = T::from(j + 1).unwrap() * n_row[j] + n_row[j - 1];
        }
    }
    n_row[k - 1]
}
