use std::fmt::{Display, Formatter};

/// An error that may occur during ligthmap generation.
#[derive(Debug)]
pub enum LightmapGenerationError {
    /// Generation was cancelled by user.
    Cancelled,
    /// An index of a vertex in a triangle is out of bounds.
    InvalidIndex,
}

impl Display for LightmapGenerationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            LightmapGenerationError::Cancelled => {
                write!(f, "Lightmap generation was cancelled by the user.")
            }
            LightmapGenerationError::InvalidIndex => {
                write!(f, "An index of a vertex in a triangle is out of bounds.")
            }
        }
    }
}
