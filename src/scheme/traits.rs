pub trait UfoParameters<const L: usize, const M: usize> {
    const A_SCHEDULE: [[usize; L]; M];
    const B_SCHEDULE: [[usize; L]; M];
}

pub trait SignatureScheme {
    type SecretKey;
    type PublicKey;
    type Message;
    type Signature;

    fn keygen(&self) -> (Self::SecretKey, Self::PublicKey);

    fn sign(&self, sk: &Self::SecretKey, m: &Self::Message) -> Self::Signature;

    fn verify(
        &self,
        pk: &Self::PublicKey,
        m: &Self::Message,
        s: &Self::Signature,
    ) -> Result<(), ()>;
}
