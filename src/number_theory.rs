use num::bigint::RandBigInt;
use num::BigInt;
use num::BigUint;
use num::Zero;
use num_traits::Pow;
use rand::distributions::uniform::SampleUniform;
use rand::Rng;

pub trait Gcd {
    ///
    /// # Examples
    ///
    /// ```
    /// use num::BigInt;
    /// use num::Integer;
    /// use num::{One, Zero};
    /// use num_theory_utils::Gcd;
    ///
    /// assert_eq!(BigInt::from(44u64), BigInt::from(2024u64).gcd_euclid(&BigInt::from(748u64)));
    /// ```

    /// Determine [greatest common divisor](https://en.wikipedia.org/wiki/Greatest_common_divisor)
    /// using the [Euclidean algorithm](https://en.wikipedia.org/wiki/Euclidean_algorithm).
    fn gcd_euclid(&self, other: &Self) -> Self;
}

impl Gcd for BigInt {
    ///
    /// GCD Calculator - The Euclidean Algorithm
    /// Input: A pair of integers a and b, not both equal to zero
    /// Output: gcd(a, b)
    ///
    fn gcd_euclid(&self, other: &BigInt) -> BigInt {
        let zero = BigInt::from(0u64);
        let mut a = self.clone();
        let mut b = other.clone();
        if a.is_zero() {
            return b;
        }

        if b.is_zero() {
            return a;
        }

        let gcd: BigInt;
        if b > a {
            gcd = b.gcd_euclid(&a);
        } else {
            let mut r: BigInt = &a % &b;
            while r > zero {
                // let q = &a / &b;
                r = &a % &b;

                if r != zero {
                    a = b;
                    b = r.clone();
                }
            }

            gcd = b;
        }

        gcd
    }
}

pub fn log_base_2(x: &BigInt) -> Result<f64, String> {
    if let Ok(log_x) = log_base_10(x) {
        let log_2 = log_base_10(&BigInt::from(2u32)).unwrap();
        let res = log_x / log_2;
        return Ok(res);
    }
    Err("Failed to calculate log base 2".into())
}

pub fn log_base_10(x: &BigInt) -> Result<f64, String> {
    use std::cmp::Ordering;
    let zero = BigInt::zero();
    let x: BigUint = match x.cmp(&zero) {
        Ordering::Less => (-x).to_biguint().unwrap(),
        Ordering::Greater => x.to_biguint().unwrap(),
        Ordering::Equal => return Err("abs_log(0)".to_string()),
    };
    let x: Vec<u8> = x.to_bytes_le();
    use num::Float;
    const BYTES: usize = 12;
    let start = if x.len() < BYTES { 0 } else { x.len() - BYTES };
    let mut n: f64 = 0.0;
    // n will accumulate the f64 value as we go.
    for i in x.iter().skip(start) {
        n = n / 256.0 + (*i as f64);
    }
    let ln_256: f64 = (256.0).ln();
    Ok(n.ln() + ln_256 * ((x.len() - 1) as f64))
}

///
/// Returns a non-negative integer a < m that satisfies a ≡ cˣ(mod m)
/// c: base
/// e: exponent
/// m: modulus
///
pub fn modular_pow(base: &BigInt, e: &BigInt, modulus: &BigInt) -> BigInt {
    // initialization
    let (zero, one, two) = (BigInt::from(0u64), BigInt::from(1u64), BigInt::from(2u64));
    let mut exp = e.clone();
    let mut a: BigInt = BigInt::from(1u64);
    let mut s: BigInt = base % modulus;

    // Converts exponent to its binary representation
    // Go through the digits from LSB to MSB in each iteration
    // if the digit == 1, a = a * s % modulus, s = s * s
    // if digit == 0, s = s * s
    while exp > zero {
        // Extract the LSB from the exp.
        if &exp & &one == one {
            a = (a * &s) % modulus;
        }

        s = (&s * &s) % modulus;

        // Division by 2 to get the next digit
        exp /= &two;
    }

    a
}

///
/// Generate a random integer in a given range
///
pub fn generate_random_int_in_range<T: Clone + PartialOrd + SampleUniform>(a: &T, b: &T) -> T {
    let mut rng = rand::thread_rng();
    // return a random BigInt between a and b
    rng.gen_range(a.clone()..b.clone())
}

///
/// Generate a random big integer in a given range
///
pub fn generate_random_bigint_in_range(a: &BigInt, b: &BigInt) -> BigInt {
    let mut rng = rand::thread_rng();
    // return a random BigInt between a and b
    rng.gen_bigint_range(a, b)
}

/// Miller-Rabin Test Step-1
/// It accepts an integer and returns a boolean value
/// 1. Express n - 1 as 2ᶠm
pub fn miller_rabin_primality(n: &BigInt) -> bool {
    let (zero, one, two) = (BigInt::from(0u64), BigInt::from(1u64), BigInt::from(2u64));
    let three = BigInt::from(3u64);
    if n <= &one || n == &BigInt::from(4u64) {
        return false;
    }
    if n <= &three {
        return true;
    }

    let mut d: BigInt = n - &one;
    // Express n - 1 as 2ᶠ.m
    while &d % 2 == zero {
        d = &d / 2;
    }
    // d = (n - 1) / 2ᶠ
    // we will iterate 5 times to ensure the result is correct
    for _ in 0..5 {
        let mut m = d.clone();
        // Randomly generate a base: a such that 1 < a < n - 1
        let a: BigInt = generate_random_int_in_range(&two, &(n - 1));

        // Calculate x ≡ a^d(mod n)
        let mut x = modular_pow(&a, &m, n);

        // if x ≡ ±1 (mod n), return true
        if x == one || x == n - 1 {
            return true;
        }

        // if x ≢ ±1 (mod n), while d != n-1 .
        // d was obtained by repeated division of (m - 1) by 2.
        // multiplying it with 2 repeatedly until it equals (m - 1)
        while m != n - 1 {
            // sqaure x - This is a^((2^j)m)(mod n)
            x = modular_pow(&x, &two, n);

            // if x ≡ -1 (mod n) the input number is probably prime
            if x == n - 1 {
                return true;
            }

            // if x ≡ -1 (mod n), then x is a factor of n
            if x == one {
                return false;
            }

            // multiplication by 2
            m *= &two;
        }
    }

    false
}

/// Miller-Rabin Test - Step 2
///
fn miller_test(d: &BigInt, n: &BigInt) -> bool {
    let (_zero, one, two) = (BigInt::from(0u64), BigInt::from(1u64), BigInt::from(2u64));
    let mut d = d.clone();
    // Randomly generate a base: a such that 1 < a < n - 1
    let a: BigInt = generate_random_int_in_range(&two, &(n - 1));

    // Calculate x ≡ a^d(mod n)
    let mut x = modular_pow(&a, &d, n);

    // if x ≡ ±1 (mod n), return true
    if x == one || x == n - 1 {
        return true;
    }

    // if x ≢ ±1 (mod n), while d != n-1 .
    // d was obtained by repeated division of (m - 1) by 2.
    // multiplying it with 2 repeatedly until it equals (m - 1)
    while d != n - 1 {
        // sqaure x - This is a^((2^j)m)(mod n)
        x = modular_pow(&x, &two, n);

        // if x ≡ -1 (mod n) the input number is probably prime
        if x == n - 1 {
            return true;
        }

        // if x ≡ -1 (mod n), then x is a factor of n
        if x == one {
            return false;
        }

        // multiplication by 2
        d *= &two;
    }

    false
}

// https://cstheory.stackexchange.com/questions/16487/find-next-prime
pub fn next_prime(n: &BigInt) -> Option<BigInt> {
    let log_n = log_base_10(n).unwrap();
    // TODO: What happens if log_n * log_n is u128?
    let upper_bound = n + BigInt::from((log_n * log_n) as u64);
    num::range(n.clone() + 1, upper_bound).find(miller_rabin_primality)
}

pub fn generate_prime_number(bit_size: usize) -> BigInt {
    let lb_exp: usize = bit_size - 1;
    let ub_exp: usize = bit_size;
    let lb = BigInt::from(2u32).pow(lb_exp);
    let ub = BigInt::from(2u32).pow(ub_exp);
    loop {
        let randp = generate_random_bigint_in_range(&lb, &ub);
        if log_base_2(&randp).unwrap().ceil() as usize == bit_size {
            if miller_rabin_primality(&randp) {
                if &randp % 4 == BigInt::from(3u32) {
                    return randp;
                }
            } else {
                let randp = next_prime(&randp).unwrap();
                if &randp % 4 == BigInt::from(3u32) {
                    return randp;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_modular_power() {
        let result = modular_pow(
            &BigInt::from(2u64),
            &BigInt::from(2u64),
            &BigInt::from(10u64),
        );
        assert_eq!(result, BigInt::from(4u64));
    }

    #[test]
    fn test_gcd_euclid() {
        let a = BigInt::from(100u64);
        let result = a.gcd_euclid(&BigInt::from(76u64));
        assert_eq!(result, BigInt::from(4u64));
        assert_eq!(
            BigInt::from(44u64),
            BigInt::from(2024u64).gcd_euclid(&BigInt::from(748u64))
        );
    }

    #[test]
    fn test_next_prime() {
        assert_eq!(BigInt::from(11u64), next_prime(&BigInt::from(10)).unwrap());
        assert_eq!(BigInt::from(13u64), next_prime(&BigInt::from(11)).unwrap());
        assert_ne!(BigInt::from(51u64), next_prime(&BigInt::from(49)).unwrap());
        assert_eq!(BigInt::from(101u64), next_prime(&BigInt::from(97)).unwrap());
        assert_eq!(
            "35241112858322930737959052235806536631839484380166598440311978859491644932877"
                .parse::<BigInt>()
                .unwrap(),
            next_prime(
                &"35241112858322930737959052235806536631839484380166598440311978859491644932747"
                    .parse::<BigInt>()
                    .unwrap()
            )
            .unwrap()
        );
    }

    #[test]
    fn test_generate_prime() {
        generate_prime_number(256);
    }
}
