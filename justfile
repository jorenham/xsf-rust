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

[working-directory("tools")]
lint-py:
    uv run ruff check
    uv run ruff format --check
    uv run pyrefly check

doc *args:
    RUSTDOCFLAGS="--html-in-header katex.html" cargo +nightly doc --all-features --no-deps {{ args }}
