default:
    @just --list

build *args:
    cargo build {{ args }}

test *args:
    cargo test {{ args }}

check *args:
    cargo check {{ args }}

clean *args:
    cargo clean {{ args }}

clippy *args:
    cargo clippy --all-targets --all-features {{ args }}

fmt *args:
    cargo fmt {{ args }}

accuracy *args:
    uv run tools/accuracy/accuracy.py {{ args }}

lint-py:
    uvx ruff@0.17.3 check tools
    uvx ruff@0.17.3 format --check tools
    uvx --with mpmath==1.4.1 basedpyright@1.40.2 -p tools

doc *args:
    RUSTDOCFLAGS="--html-in-header katex.html" cargo +nightly doc --all-features --no-deps {{ args }}
