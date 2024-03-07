use std::fmt::{Display, Formatter};

/// An error that may occur during ligth map generation.
#[derive(Debug)]
pub enum Error {
    /// An index of a vertex in a triangle is out of bounds.
    InvalidIndex,
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidIndex => {
                write!(f, "An index of a vertex in a triangle is out of bounds.")
            }
        }
    }
}
