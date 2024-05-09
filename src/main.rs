use std::fmt::Debug;

use anyhow::Result;
use num::BigInt;

use crate::{
    algebra::PrimeFieldConfig,
    ecdisplay::EcDisplay,
    elgamal::{DiffieKeyExchange, EncryptionScheme},
    elliptic_curves::{CurveParams, CurvePoint, EllipticCurve},
    number_theory::generate_prime_number,
    utf8::{bytes_to_utf8, utf8_to_bytes},
};

mod algebra;
mod ecdisplay;
mod elgamal;
mod elliptic_curves;
mod number_theory;
mod utf8;

fn main() {
    // Generation of a 256-bit curve
    if let Err(e) = test_256bit_ecc() {
        eprintln!("{}", e)
    }

    // Generation of a 521-bit curve
    if let Err(e) = test_521bit_ecc() {
        eprintln!("{}", e)
    }
}

fn test_256bit_ecc() -> Result<()> {
    // Step1: Generate a 256-bit prime number
    let prime = generate_prime_number(256);
    println!("\nPrime Number: {}", &prime);
    // Step2: Create a prime field which will be inherited
    // in further operations
    let q = prime_field!(prime.clone());
    // Step3: Generate the curve parameters
    let curve1_params = CurveParams::init(&q);
    // Step4: Initialise the curve with the generator point, G
    // using the curve parameters. ``None", ``None" below
    // will force generation of new values for G.
    let curve1 = CurvePoint::new(None, None, curve1_params);
    println!("Generated Curve:");
    println!("{}", curve1.format_curve().unwrap());

    // Step5: Generate Private and Public Keys for Alice
    println!("Generating Private Key and Public Key for Alice:");
    let alice = DiffieKeyExchange::new(curve1.clone());
    let alice_private_key = alice
        .get_private_key()
        .expect("Retrieving private key failed!");
    println!("  Alice's Private Key: \n\t{}", alice_private_key.clone().0);
    let alice_public_key = alice
        .get_public_key()
        .expect("Retrieving public key failed");
    println!("  Alice's Public Key: \n\t{}", alice_public_key);

    // Step6: Generate Private and Public Keys for Bob
    println!("\nGenerating Private Key and Public Key for Bob:");
    let bob = DiffieKeyExchange::new(curve1.clone());
    let bob_private_key = bob
        .get_private_key()
        .expect("Retrieving private key failed");
    println!("  Bob's Private Key: \n\t{}", bob_private_key.clone().0);
    let bob_public_key = bob.get_public_key().expect("Retrieving public key failed");
    println!("  Bob's Public Key: \n\t{}\n", bob_public_key);

    // let msg1 = "Hello World";
    let msg1 = "London is vibrant city!!!! River.";
    println!("Original message at Alice's: \n\t{}\n", msg1);
    // Step7: @Alice -Convert the utf8 input string to hexadecimal bytes
    let msg1 = utf8_to_bytes(msg1).unwrap();
    // Step8: Encode the message into EC point using Koblitz algm
    let alice_encoded_message1 = alice.koblitz_point_encode(&msg1)?;
    println!(
        "Koblitz encoded message at Alice's: \n{}\n",
        alice_encoded_message1.clone()
    );

    // Step9: @Alice - Encrypt the message using Bob's public key and Alice's
    // private key
    let alice_encrypt_msg1 = alice_encoded_message1 + alice_private_key.0 * bob_public_key.0;
    println!("Encrypted message at Alice's: \n{}\n", alice_encrypt_msg1);

    // Step10: @Bob - Decrypt ciphertext from Alice using Bob's private key
    // Alice's public key
    let bob_decrypt_msg1 = alice_encrypt_msg1 + bob_private_key.0 * -alice_public_key.0;
    println!("Decrypted message at Bob's: \n{}\n", bob_decrypt_msg1);

    // Step11: @Bob - Decode the decrypted message
    let bob_decoded_message1 = bob.koblitz_point_decode(&bob_decrypt_msg1).unwrap();
    println!(
        "Koblitz decoded message at Bob's: \n{:?}\n",
        bob_decoded_message1.clone()
    );

    // Step12: @Bob - convert the hexademical bytes into utf8 string
    let recovered_msg = bytes_to_utf8(bob_decoded_message1)?;
    println!("Recovered message at Bob's: \n\t{}\n", recovered_msg);

    Ok(())
}

fn test_521bit_ecc() -> Result<()> {
    let prime = generate_prime_number(521);
    println!("\nPrime Number: {}", &prime);
    let q = prime_field!(prime.clone());
    let curve1_params = CurveParams::init(&q);
    let curve1 = CurvePoint::new(None, None, curve1_params);
    println!("Generated Curve:");
    println!("{}", curve1.format_curve().unwrap());

    println!("Generating Private Key and Public Key for Alice:");
    let alice = DiffieKeyExchange::new(curve1.clone());
    let alice_private_key = alice
        .get_private_key()
        .expect("Retrieving private key failed!");
    println!("  Alice's Private Key: \n\t{}", alice_private_key.clone().0);
    let alice_public_key = alice
        .get_public_key()
        .expect("Retrieving public key failed");
    println!("  Alice's Public Key: \n\t{}", alice_public_key);

    println!("\nGenerating Private Key and Public Key for Bob:");
    let bob = DiffieKeyExchange::new(curve1.clone());
    let bob_private_key = bob
        .get_private_key()
        .expect("Retrieving private key failed");
    println!("  Bob's Private Key: \n\t{}", bob_private_key.clone().0);
    let bob_public_key = bob.get_public_key().expect("Retrieving public key failed");
    println!("  Bob's Public Key: \n\t{}\n", bob_public_key);

    // let msg1 = "Hello World";
    let msg1 = "Sherlock Holmes of Scotland Yards";
    println!("Original message at Alice's: \n\t{}\n", msg1);
    let msg1 = utf8_to_bytes(msg1).unwrap();
    let alice_encoded_message1 = alice.koblitz_point_encode(&msg1)?;
    println!(
        "Koblitz encoded message at Alice's: \n{}\n",
        alice_encoded_message1.clone()
    );

    let alice_encrypt_msg1 = alice_encoded_message1 + alice_private_key.0 * bob_public_key.0;
    println!("Encrypted message at Alice's: \n{}\n", alice_encrypt_msg1);

    let bob_decrypt_msg1 = alice_encrypt_msg1 + bob_private_key.0 * -alice_public_key.0;
    println!("Decrypted message at Bob's: \n{}\n", bob_decrypt_msg1);

    let bob_decoded_message1 = bob.koblitz_point_decode(&bob_decrypt_msg1).unwrap();
    println!(
        "Koblitz decoded message at Bob's: \n{:?}\n",
        bob_decoded_message1.clone()
    );

    let recovered_msg = bytes_to_utf8(bob_decoded_message1)?;
    println!("Recovered message at Bob's: \n\t{}\n", recovered_msg);

    Ok(())
}
