/// `$strategy`, or else `arg()`
macro_rules! strategy {
    () => {
        $crate::properties::strategies::arg()
    };
    ($strategy:expr) => {
        $strategy
    };
}

/// A test `$name` that passes `$call` to `$check`, unless `$assume` is false, e.g.
/// `property!(no_panic, beta(a, b in positive()))` checks `no_panic(xsf::beta(a, b))` with `a` from
/// `arg()`
macro_rules! property {
    (
        $check:path,
        $name:ident($($arg:ident $(in $strategy:expr)?),+) $(if $assume:expr)? => $call:expr
    ) => {
        proptest::proptest! {
            #![proptest_config($crate::properties::config())]
            #[test]
            fn $name($($arg in strategy!($($strategy)?)),+) {
                $(proptest::prop_assume!($assume);)?
                $check($call)?;
            }
        }
    };
    ($check:path, $name:ident($($arg:ident $(in $strategy:expr)?),+) $(if $assume:expr)?) => {
        property!(
            $check, $name($($arg $(in $strategy)?),+) $(if $assume)? => xsf::$name($($arg),+)
        );
    };
}
