use crate::params::{Level1Config, Level3Config, Level5Config};
use crate::scheme::UFOs;

pub mod field;
pub mod matrix;
pub mod params;
pub mod random;
pub mod scheme;

type UfosLevel1 = UFOs<31, 61, Level1Config, 5, 5, 1, 1, 4>;

type UfosLevel3 = UFOs<61, 89, Level3Config, 5, 5, 1, 1, 4>;

type UfosLevel5 = UFOs<61, 127, Level5Config, 5, 5, 1, 1, 4>;
