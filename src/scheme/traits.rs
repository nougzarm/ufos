use crate::field::{
    traits::ExtensionConfig,
    prime_field::Fp,
};

pub trait UFOs {
    const P: u8;
    const E: usize;

}