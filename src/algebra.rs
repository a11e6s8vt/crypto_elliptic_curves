use num::{range, BigInt, One, Zero};
use num_traits::Pow;
use std::{
    fmt::Debug,
    marker::PhantomData,
    ops::{Add, Div, DivAssign, Mul, Neg, Rem, Sub},
};

use crate::number_theory::{generate_random_int_in_range, modular_pow};

// Presents a field.
pub trait PrimeFieldConfig<F>
where
    F: Clone
        + 'static
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Neg<Output = F>
        + Sub<F, Output = F>
        + Div<F, Output = F>
        + Rem<Output = F>,
{
    fn p(&self) -> F;
}

#[macro_export]
macro_rules! prime_field {
    ($N:expr) => {{
        // define the struct
        #[derive(Debug, Clone)]
        pub struct PrimeField(BigInt);

        impl PrimeFieldConfig<BigInt> for PrimeField {
            fn p(&self) -> BigInt {
                self.0.clone()
            }
        }

        // construct an instance of the struct and return it
        PrimeField($N)
    }};
}

// Represents the field Z_p for the given prime number.
#[derive(Clone, Debug)]
pub struct FieldElement<'a, P>(pub BigInt, pub &'a P, std::marker::PhantomData<P>)
where
    P: PrimeFieldConfig<BigInt> + Clone + Debug;

// impl<'a, P: PrimeFieldConfig<BigInt>> Clone for FieldElement<'a, P> {
//     fn clone(&self) -> Self {}
// }

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> FieldElement<'a, P> {
    pub fn new(fp: &'a P, n: Option<&BigInt>) -> FieldElement<'a, P> {
        let one = BigInt::one();
        let n = n.unwrap_or(&one);
        let result = n % fp.p();
        if result < BigInt::zero() {
            FieldElement(result + fp.p(), fp, PhantomData)
        } else {
            FieldElement(result, fp, PhantomData)
        }
    }

    ///
    /// The Extended Euclidean Algorithm to find modular inverse
    ///
    pub fn inv(self) -> Option<FieldElement<'a, P>> {
        if self.0 == BigInt::one() {
            return Some(self);
        }

        let mut u = vec![self.1.p(), BigInt::one(), BigInt::zero()];
        let mut v = vec![self.0, BigInt::zero(), BigInt::one()];

        while v[0] > BigInt::zero() {
            let q = u[0].clone() / v[0].clone();
            let i = u.get(0).unwrap().clone() - q.clone() * v.get(0).unwrap().clone();
            let j = u.get(1).unwrap().clone() - q.clone() * v.get(1).unwrap().clone();
            let k = u.get(2).unwrap().clone() - q.clone() * v.get(2).unwrap().clone();
            let _ = std::mem::replace(&mut u, v.clone());
            let _ = std::mem::replace(&mut v, vec![i, j, k]);
        }

        if u[0] != BigInt::one() {
            return None;
        }

        if u[2] > BigInt::zero() {
            Some(FieldElement(u.get(2).unwrap().clone(), self.1, PhantomData))
        } else {
            Some(FieldElement(
                self.1.p() + u.get(2).unwrap().clone(),
                self.1,
                PhantomData,
            ))
        }
    }

    /// SQRT using Tonelli-Shanks Algorithm
    pub fn sqrt(&self) -> Option<FieldElement<'a, P>> {
        let (zero, one, two) = (BigInt::from(0u64), BigInt::from(1u64), BigInt::from(2u64));
        if modular_pow(&self.0, &((self.1.p() - &one) / &two), &self.1.p()) == BigInt::from(-1i32) {
            return None;
        }

        let p = self.1.p();
        if &p % &BigInt::from(4u32) == BigInt::from(3u32) {
            return Some(FieldElement::new(
                self.1,
                Some(
                    &(modular_pow(
                        &self.0,
                        &((p.clone() + one.clone()) / BigInt::from(4u32)),
                        &p,
                    )),
                ),
            ));
        } else {
            let p_minus_one = p.clone() - BigInt::one();
            let mut m = p_minus_one.clone();

            let mut s = 0;
            while &m % &two == zero {
                m = &m / &two;
                s += 1;
            }
            // finding a quadratic non-residue
            let mut c: BigInt = zero.clone();
            //loop {

            for a in range(two.clone(), p_minus_one.clone()) {
                // let a = generate_random_int_in_range(&BigInt::from(2u32), &p_minus_one);
                let g = modular_pow(&a, &(p_minus_one.clone() / two.clone()), &p);

                if g == p_minus_one {
                    c = modular_pow(&a, &m, &p);
                    break;
                }
            }

            let mut u = modular_pow(&self.0, &m, &p);
            let mut v = modular_pow(&self.0, &((&m + &one) / &two), &p);

            while u != one {
                let mut k: u32 = 0;
                let mut l: u32 = 0;

                loop {
                    let c_raised_2_k = modular_pow(&c.clone(), &BigInt::from(2u32).pow(k), &p);
                    if FieldElement::new(self.1, Some(&c_raised_2_k))
                        == FieldElement::new(self.1, Some(&p_minus_one))
                    {
                        break;
                    } else {
                        k += 1;
                    }
                }

                loop {
                    let u_raised_2_l = modular_pow(&u, &BigInt::from(2u32).pow(l), &p);
                    if u_raised_2_l != p_minus_one {
                        l += 1;
                    } else {
                        break;
                    }
                }

                let r_k = l + 1;
                u = FieldElement::new(
                    self.1,
                    Some(&(u.clone() * modular_pow(&c, &BigInt::from(2u32).pow(s - r_k), &p))),
                )
                .0;

                v = FieldElement::new(
                    self.1,
                    Some(&(v.clone() * modular_pow(&c, &BigInt::from(2u32).pow(s - r_k - 1), &p))),
                )
                .0;
            }

            Some(FieldElement::new(self.1, Some(&v)))
        }
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> PartialEq<FieldElement<'a, P>>
    for FieldElement<'a, P>
{
    fn eq(&self, other: &FieldElement<'a, P>) -> bool {
        if self.1.p() == other.1.p() && self.0 == other.0 {
            return true;
        }

        false
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Add<FieldElement<'a, P>>
    for FieldElement<'a, P>
{
    type Output = FieldElement<'a, P>;

    fn add(self, rhs: FieldElement<P>) -> Self::Output {
        Self::new(self.1, Some(&(self.0 + rhs.0)))
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Mul<FieldElement<'a, P>>
    for FieldElement<'a, P>
{
    type Output = FieldElement<'a, P>;

    fn mul(self, rhs: FieldElement<'a, P>) -> Self::Output {
        Self::new(self.1, Some(&(self.0 * rhs.0)))
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Div<FieldElement<'a, P>>
    for FieldElement<'a, P>
{
    type Output = Option<FieldElement<'a, P>>;

    fn div(self, rhs: FieldElement<'a, P>) -> Self::Output {
        if let Some(rhs_inv) = rhs.inv() {
            let res = self * rhs_inv;
            return Some(res);
        } else {
            None
        }
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Neg for FieldElement<'a, P> {
    type Output = FieldElement<'a, P>;

    fn neg(self) -> Self::Output {
        Self::new(self.1, Some(&-self.0))
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Sub<FieldElement<'a, P>>
    for FieldElement<'a, P>
{
    type Output = FieldElement<'a, P>;

    fn sub(self, rhs: FieldElement<'a, P>) -> Self::Output {
        Self::new(self.1, Some(&(self.0 - rhs.0)))
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Pow<u32> for FieldElement<'a, P> {
    type Output = FieldElement<'a, P>;

    fn pow(self, rhs: u32) -> Self::Output {
        FieldElement::new(self.1, Some(&self.0.pow(rhs)))
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Rem for FieldElement<'a, P> {
    type Output = FieldElement<'a, P>;

    fn rem(self, rhs: Self) -> Self::Output {
        let result = self.0 % rhs.clone().0;
        if result < BigInt::zero() {
            FieldElement::new(&self.1, Some(&(result + rhs.0)))
        } else {
            FieldElement::new(&self.1, Some(&result))
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use num::BigInt;

    #[test]
    fn test_modular_add_case1() {
        let a: BigInt = BigInt::from(19u64);
        let field = prime_field!(a);

        let p = FieldElement::new(&field, Some(&BigInt::from(17u64)));
        let q = FieldElement::new(&field, Some(&BigInt::from(5u64)));
        let res = FieldElement::new(&field, Some(&BigInt::from(3u64)));
        assert_eq!(res, p + q);
    }

    #[test]
    fn test_modular_sub_case1() {
        let a: BigInt = BigInt::from(19u64);
        let field = prime_field!(a);

        let p = FieldElement::new(&field, Some(&BigInt::from(23u64)));
        let q = FieldElement::new(&field, Some(&BigInt::from(17u64)));
        let res = FieldElement::new(&field, Some(&BigInt::from(6u64)));
        assert_eq!(res, p - q);
    }

    #[test]
    fn test_modular_mul_case1() {
        let a: BigInt = BigInt::from(19u64);
        let field = prime_field!(a);

        let p = FieldElement::new(&field, Some(&BigInt::from(23u64)));
        let q = FieldElement::new(&field, Some(&BigInt::from(17u64)));
        let res = FieldElement::new(&field, Some(&BigInt::from(11u64)));
        assert_eq!(res, p * q);
    }

    #[test]
    fn test_modular_inv() {
        let a = BigInt::from(148u64);
        let field = prime_field!(a);

        let p = FieldElement::new(&field, Some(&BigInt::from(75u64)));
        let res = p.inv().unwrap();
        assert_eq!(res, FieldElement::new(&field, Some(&BigInt::from(75u64))));
    }

    #[test]
    fn test_modular_sqrt_case1() {
        let a = BigInt::from(769u64);
        let field = prime_field!(a);

        let p = FieldElement::new(&field, Some(&BigInt::from(6u64)));
        let res = p.sqrt().unwrap();
        assert_eq!(res, FieldElement::new(&field, Some(&BigInt::from(542u64))));
    }

    #[test]
    fn test_modular_sqrt_case2() {
        let a = BigInt::from(400009u64);
        let field = prime_field!(a);

        let p = FieldElement::new(&field, Some(&BigInt::from(2u64)));
        let res = p.sqrt().unwrap();
        assert_eq!(
            res,
            FieldElement::new(&field, Some(&BigInt::from(282720u64)))
        );
    }
}
