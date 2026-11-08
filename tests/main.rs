//! Tests that only use the public API, mirroring the layout of `src/`.

#![warn(clippy::pedantic)]
#![allow(
    clippy::doc_markdown,
    clippy::many_single_char_names,
    clippy::similar_names,
    // casts between the test indices, the reference data, and the argument types
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]

mod numpy;
mod properties;
mod scipy_special;
mod xsf;
