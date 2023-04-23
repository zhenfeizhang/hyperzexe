//! Main module for the HyperPlonk SNARK.

pub mod backend;
pub mod frontend;
pub use halo2_curves;

#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    InvalidSumcheck(String),
    InvalidPcsParam(String),
    InvalidPcsOpen(String),
    InvalidSnark(String),
    Serialization(String),
    Transcript(std::io::ErrorKind, String),
    NotImplemented(String),
}
