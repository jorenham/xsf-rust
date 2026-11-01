/// Complete elliptic integral of the first kind $K(m)$
///
/// This function is defined as [^DLMF]:
///
/// $$
/// K(m) = \int_0^{\pi/2} {\dd \theta \over \sqrt{1 - m \sin^2 \theta}} .
/// $$
///
/// Corresponds to [`scipy.special.ellipk`][ellipk] in SciPy
///
/// [ellipk]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipk.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/19.2.E8>, with $m = k^2$.
///
/// # See also
/// - [`ellipe`]: Complete elliptic integral of the second kind
/// - [`ellipkinc`]: Incomplete elliptic integral of the first kind
/// - [`ellipkm1`]: Complete elliptic integral of the first kind around $m = 1$
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_k")]
#[must_use]
#[inline]
pub fn ellipk(m: f64) -> f64 {
    crate::xsf::cephes::ellpk(1.0 - m)
}

/// Complete elliptic integral of the first kind around $m = 1$
///
/// This computes $K(1 - p)$, where $K$ is [`ellipk`].
///
/// Corresponds to [`scipy.special.ellipkm1`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkm1.html
///
/// # See also
/// - [`ellipk`]: Complete elliptic integral of the first kind
/// - [`ellipkinc`]: Incomplete elliptic integral of the first kind
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_k_m1")]
#[must_use]
#[inline]
pub fn ellipkm1(p: f64) -> f64 {
    crate::xsf::cephes::ellpk(p)
}

/// Incomplete elliptic integral of the first kind $F(\phi \vert m)$
///
/// This function is defined as [^DLMF]:
///
/// $$
/// F(\phi \vert m) = \int_0^\phi {\dd \theta \over \sqrt{1 - m \sin^2 \theta}} .
/// $$
///
/// Corresponds to [`scipy.special.ellipkinc`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkinc.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/19.2.E4>, with $m = k^2$.
///
/// # See also
/// - [`ellipk`]: Complete elliptic integral of the first kind
/// - [`ellipkm1`]: Complete elliptic integral of the first kind around $m = 1$
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_k_inc")]
#[must_use]
#[inline]
pub fn ellipkinc(phi: f64, m: f64) -> f64 {
    crate::xsf::cephes::ellik(phi, m)
}

/// Complete elliptic integral of the second kind $E(m)$
///
/// This function is defined as [^DLMF]:
///
/// $$
/// E(m) = \int_0^{\pi/2} \sqrt{1 - m \sin^2 \theta} \dd \theta .
/// $$
///
/// Corresponds to [`scipy.special.ellipe`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/19.2.E8>, with $m = k^2$.
///
/// # See also
/// - [`ellipk`]: Complete elliptic integral of the first kind
/// - [`ellipeinc`]: Incomplete elliptic integral of the second kind
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_e")]
#[must_use]
#[inline]
pub fn ellipe(m: f64) -> f64 {
    crate::xsf::cephes::ellpe(m)
}

/// Incomplete elliptic integral of the second kind $E(\phi \vert m)$
///
/// This function is defined as [^DLMF]:
///
/// $$
/// E(\phi \vert m) = \int_0^\phi \sqrt{1 - m \sin^2 \theta} \dd \theta .
/// $$
///
/// Corresponds to [`scipy.special.ellipeinc`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipeinc.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/19.2.E5>, with $m = k^2$.
///
/// # See also
/// - [`ellipe`]: Complete elliptic integral of the second kind
/// - [`ellipk`]: Complete elliptic integral of the first kind
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_e_inc")]
#[must_use]
#[inline]
pub fn ellipeinc(phi: f64, m: f64) -> f64 {
    crate::xsf::cephes::ellie(phi, m)
}

/// Jacobi elliptic functions
///
/// Corresponds to [`scipy.special.ellipj`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipj.html
///
/// # Arguments
/// - `u` - Real argument
/// - `m` - Parameter, `m <= 1`
///
/// # Returns
/// - *sn(u | m)* - sine amplitude
/// - *cn(u | m)* - cosine amplitude
/// - *dn(u | m)* - delta amplitude
/// - phase *φ* s.t. *sn(u | m) = sin(φ)* and *cn(u | m) = cos(φ)*
///
/// # See also
/// - [`ellipk`] - Complete elliptic integral of the first kind
/// - [`ellipe`] - Complete elliptic integral of the second kind
/// - [`ellipkinc`] - Incomplete elliptic integral of the first kind
/// - [`ellipeinc`] - Incomplete elliptic integral of the second kind
#[doc(alias = "ellip_j")]
#[must_use]
#[inline]
pub fn ellipj(u: f64, m: f64) -> (f64, f64, f64, f64) {
    crate::xsf::cephes::ellpj(u, m)
}
