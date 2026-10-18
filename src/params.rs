use crate::field::{prime_field::Fp, traits::ExtensionConfig};

/// Level I of UFOs : P = 31, E = 61
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct UfoLevel1Config;

impl ExtensionConfig<Fp<31>, 61> for UfoLevel1Config {
    const REDUCE_COEFFS: [Fp<31>; 61] = {
        let mut coeffs = [Fp::new(0); 61];
        coeffs[0] = Fp::new(1); 
        coeffs[1] = Fp::new(1); 
        coeffs
    };
}

pub type F31_61 = crate::field::extension::Fq<31, 61, UfoLevel1Config>;
