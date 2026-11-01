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
    uv run --project tools --locked --only-group lint ruff check tools
    uv run --project tools --locked --only-group lint ruff format --check tools
    uv run --project tools --locked --only-group lint basedpyright -p tools

doc *args:
    RUSTDOCFLAGS="--html-in-header katex.html" cargo +nightly doc --all-features --no-deps {{ args }}
