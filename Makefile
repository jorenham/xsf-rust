doc:
	RUSTDOCFLAGS="--html-in-header katex.html" cargo +nightly doc --all-features --no-deps
