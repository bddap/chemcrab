default:
    @just --list

check:
    cargo fmt --check
    cargo clippy -- -D warnings
    RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
    cargo test
