doc:
	RUSTDOCFLAGS="--html-in-header katex.html" cargo +nightly doc --all-features --no-deps
	echo "./target/docs/math-in-docs/index.html"
