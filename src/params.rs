use crate::field::{prime_field::Fp, traits::ExtensionConfig};
use crate::scheme::traits::UfoParameters;

/// Level I of UFOs : P = 31, E = 61
#[derive(Copy, Clone, Default, PartialEq, Debug)]
pub struct Level1Config;

impl ExtensionConfig<61> for Level1Config {
    type Field = Fp<31>;

    const REDUCE_COEFFS: [Self::Field; 61] = {
        let mut coeffs = [Self::Field::new(0); 61];
        coeffs[0] = Fp::from_canonical(3);
        coeffs[1] = Fp::from_canonical(25);
        coeffs[2] = Fp::from_canonical(1);
        coeffs
    };
}

impl UfoParameters<5, 1> for Level1Config {
    const A_SCHEDULE: [[usize; 5]; 1] = [[0, 7, 13, 15, 35]];
    const B_SCHEDULE: [[usize; 5]; 1] = [[9, 54, 29, 37, 47]];
}

/// Level III of UFOs : P = 61, E = 89
#[derive(Copy, Clone, Default, PartialEq, Debug)]
pub struct Level3Config;

impl ExtensionConfig<89> for Level3Config {
    type Field = Fp<61>;

    const REDUCE_COEFFS: [Self::Field; 89] = {
        let mut coeffs = [Self::Field::new(0); 89];
        coeffs[0] = Fp::from_canonical(2);
        coeffs[1] = Fp::from_canonical(30);
        coeffs
    };
}

impl UfoParameters<5, 1> for Level3Config {
    const A_SCHEDULE: [[usize; 5]; 1] = [[0, 28, 40, 52, 53]];
    const B_SCHEDULE: [[usize; 5]; 1] = [[13, 43, 61, 84, 76]];
}

/// Level V of UFOs : P = 61, E = 127
#[derive(Copy, Clone, Default, PartialEq, Debug)]
pub struct Level5Config;

impl ExtensionConfig<127> for Level5Config {
    type Field = Fp<61>;

    const REDUCE_COEFFS: [Self::Field; 127] = {
        let mut coeffs = [Self::Field::new(0); 127];
        coeffs[0] = Fp::from_canonical(46);
        coeffs[6] = Fp::from_canonical(60);
        coeffs
    };
}

impl UfoParameters<5, 1> for Level5Config {
    const A_SCHEDULE: [[usize; 5]; 1] = [[0, 22, 29, 53, 94]];
    const B_SCHEDULE: [[usize; 5]; 1] = [[2, 83, 74, 79, 111]];
}

pub type F31_61 = crate::field::extension::Fq<31, 61, Level1Config>;
pub type F61_89 = crate::field::extension::Fq<61, 89, Level3Config>;
pub type F61_127 = crate::field::extension::Fq<61, 127, Level5Config>;
