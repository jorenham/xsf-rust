//! No panics, e.g. from the overflow checks of debug builds, over the full `f64` range

use num_complex::c64;
use proptest::prelude::*;

#[allow(clippy::unnecessary_wraps)] // a `property!` check
fn no_panic<T>(value: T) -> Result<(), TestCaseError> {
    std::hint::black_box(value);
    Ok(())
}

property!(no_panic, hyp2f1(a, b, c, x));
property!(no_panic, hyp2f1_complex(a, b, c, re, im) => xsf::hyp2f1(a, b, c, c64(re, im)));
property!(no_panic, hyp1f1(a, b, x));
property!(no_panic, hypu(a, b, x));
property!(no_panic, hyp0f1(b, x));

property!(no_panic, elliprc(x, y));
property!(no_panic, elliprd(x, y, z));
property!(no_panic, elliprf(x, y, z));
property!(no_panic, elliprg(x, y, z));
property!(no_panic, elliprj(x, y, z, p));

property!(no_panic, betainc(a, b, x));
property!(no_panic, betaincinv(a, b, y));
property!(no_panic, betaincc(a, b, x));
property!(no_panic, betainccinv(a, b, y));
property!(no_panic, gammainc(a, x));
property!(no_panic, gammaincc(a, x));
property!(no_panic, gammaincinv(a, p));
property!(no_panic, gammainccinv(a, p));

property!(no_panic, fdtr(a, b, x));
property!(no_panic, stdtr(nu, x));
property!(no_panic, stdtri(nu, p));
property!(no_panic, bdtr(k, n in any::<i32>(), p));
property!(no_panic, nbdtr(k in any::<i32>(), n in any::<i32>(), p));
property!(no_panic, pdtr(k, m));
property!(no_panic, smirnov(n in any::<i32>(), x));
property!(no_panic, owens_t(h, a));

property!(no_panic, bessel_j(v, x));
property!(no_panic, bessel_y(v, x));
property!(no_panic, bessel_i(v, x));
property!(no_panic, bessel_k(v, x));
property!(no_panic, bessel_j_complex(v, re, im) => xsf::bessel_j(v, c64(re, im)));
property!(no_panic, bessel_k_complex(v, re, im) => xsf::bessel_k(v, c64(re, im)));
property!(no_panic, struve_h(v, x));
property!(no_panic, struve_l(v, x));

property!(no_panic, pbdv(v, x));
property!(no_panic, pbvv(v, x));
property!(no_panic, pbwa(a, x));
property!(no_panic, wright_bessel(a, b, x));

property!(no_panic, mathieu_a(m in 0..10_000_u32, q));
property!(no_panic, mathieu_b(m in 0..10_000_u32, q));
property!(no_panic, mathieu_cem(m in 0..10_000_u32, q, x));
property!(no_panic, mathieu_sem(m in 0..10_000_u32, q, x));

property!(no_panic, hahn_q(n in 0..1000_u32, big_n in any::<u32>(), alpha, beta, x));
