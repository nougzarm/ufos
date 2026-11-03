use rand::RngCore;
use sha3::digest::{ExtendableOutput, Update, XofReader};

pub fn sample_secret_seed() -> [u8; 32] {
    let mut seed = [0u8; 32];
    rand::rngs::OsRng.fill_bytes(&mut seed);
    seed
}

pub fn sample_public_seed(secret_seed: &[u8; 32]) -> [u8; 16] {
    let mut hasher = sha3::Shake256::default();
    hasher.update(b"UFOs_Public_Seed");
    hasher.update(secret_seed);

    let mut reader = hasher.finalize_xof();
    let mut public_seed = [0u8; 16];
    reader.read(&mut public_seed);

    public_seed
}

pub fn sample_salt() -> [u8; 32] {
    let mut salt = [0u8; 32];
    rand::rngs::OsRng.fill_bytes(&mut salt);
    salt
}
