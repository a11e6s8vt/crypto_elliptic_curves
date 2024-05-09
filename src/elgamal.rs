use std::fmt::Debug;
use std::fmt::Display;

use anyhow::anyhow;
use anyhow::Result;
use num::{range, BigInt, Zero};
use num_traits::FromBytes;
use num_traits::Pow;

use crate::{
    algebra::{FieldElement, PrimeFieldConfig},
    elliptic_curves::{CurvePoint, EllipticCurve},
    number_theory::generate_random_int_in_range,
};

/// An ECDH private key
#[derive(Clone, Debug)]
pub struct PrivateKey(pub BigInt);

/// An ECDH public key.
#[derive(Debug, Clone, PartialEq)]
pub struct PublicKey<'a, P>(pub CurvePoint<'a, P>)
where
    P: PrimeFieldConfig<BigInt> + Clone + Debug;

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Display for PublicKey<'a, P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "    {}", self.0)
    }
}
/// A trait to describe the types, methods and functions of a key-exhange for a curve
pub trait EncryptionScheme<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> {
    fn new(curve: CurvePoint<'a, P>) -> Self;
    // Encode a plaintext message
    fn koblitz_point_encode(&'a self, message: &[u8]) -> anyhow::Result<CurvePoint<'a, P>>;
    // Decode a plaintext message
    fn koblitz_point_decode(&self, msg_ec_point: &CurvePoint<'a, P>) -> Result<Vec<u8>>;
}

/// A struct that represents the ECDH implementation for the p-256 curve
#[derive(Clone, Debug)]
pub struct DiffieKeyExchange<'a, P>
where
    P: PrimeFieldConfig<BigInt> + Clone + Debug,
{
    private_key: Option<PrivateKey>,
    public_key: Option<PublicKey<'a, P>>,
    curve: CurvePoint<'a, P>,
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> DiffieKeyExchange<'a, P> {
    pub fn get_private_key(&self) -> Option<PrivateKey> {
        self.private_key.clone()
    }

    pub fn get_public_key(&self) -> Option<PublicKey<'a, P>> {
        self.public_key.clone()
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> EncryptionScheme<'a, P>
    for DiffieKeyExchange<'a, P>
{
    fn new(curve: CurvePoint<'a, P>) -> Self {
        let underlying_field = curve.underlying_field();
        let p = underlying_field.p();

        // Private Key calculation
        let sqrt_p = p.sqrt();
        let lower_bound = sqrt_p.clone();
        let upper_bound = p.clone() - sqrt_p.clone();
        let private_key = generate_random_int_in_range(&(lower_bound + 1), &upper_bound);
        let private_key = PrivateKey(private_key);

        // Public Key computation
        let pubkey_point = private_key.clone().0 * curve.clone();
        let public_key = PublicKey(pubkey_point);

        Self {
            private_key: Some(private_key),
            public_key: Some(public_key),
            curve,
        }
    }

    fn koblitz_point_encode(&'a self, message: &[u8]) -> anyhow::Result<CurvePoint<'a, P>> {
        let underlying_field = self.curve.underlying_field();
        let p = underlying_field.p();
        let msg_lower_bound = BigInt::zero();
        let msg_upper_bound = (p.clone() / 1000) - 1;

        let message: BigInt = FromBytes::from_be_bytes(message);

        // Below code will make the program fail if the condition do not satisfies
        // TODO: Split the message into multiple chunks and then encrypt
        if !(message > msg_lower_bound && message < msg_upper_bound) {
            return Err(anyhow!(
                "[ERR] - Error encoding message! Long message. I can only encrypt smaller message!!"
            ));
        }

        for i in range(BigInt::zero(), BigInt::from(1000u32)) {
            let x = (BigInt::from(1000u32) * &message) + i;
            let x = FieldElement::new(underlying_field, Some(&x));
            let y_square = x.clone().pow(3u32)
                + self.curve.get_curve_coeff_a().unwrap() * x.clone()
                + self.curve.get_curve_coeff_b().unwrap().clone();
            if let Some(y) = y_square.sqrt() {
                return Ok(CurvePoint::new(
                    Some(x),
                    Some(y),
                    self.curve.get_curve_params().unwrap(),
                ));
            }
        }
        Err(anyhow!("Failed to encode the message to points on curve!!"))
    }

    fn koblitz_point_decode(&self, msg_ec_point: &CurvePoint<'a, P>) -> Result<Vec<u8>> {
        let x = msg_ec_point.get_coordinate_x().unwrap();
        let msg = x.0 / BigInt::from(1000u32);
        Ok(msg.to_bytes_be().1)
    }
}
