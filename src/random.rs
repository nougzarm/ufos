use rand::RngCore;

pub fn sample_secret_seed() -> [u8; 32] {
    let mut seed = [0u8; 32];
    rand::rngs::OsRng.fill_bytes(&mut seed);
    seed
}

pub fn sample_public_seed() -> [u8; 16] {
    let mut seed = [0u8; 16];
    rand::rngs::OsRng.fill_bytes(&mut seed);
    seed
}

pub fn sample_salt() -> [u8; 32] {
    let mut salt = [0u8; 32];
    rand::rngs::OsRng.fill_bytes(&mut salt);
    salt
}
