use std::fmt::{Debug, Display};
use std::ops::{Add, Mul, Neg};

use crate::algebra::{FieldElement, PrimeFieldConfig};
use crate::ecdisplay::EcDisplay;
use crate::number_theory::{generate_random_int_in_range, modular_pow};
use crate::prime_field;
use fmtastic::Superscript;
use num::{BigInt, One, Zero};
use num_traits::Pow;

#[derive(Clone, Debug)]
pub struct CurveParams<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> {
    coeff_a: FieldElement<'a, P>,
    coeff_b: FieldElement<'a, P>,
    disc: FieldElement<'a, P>,
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> CurveParams<'a, P> {
    pub fn init(field: &'a P) -> Self {
        fn calc_curve_params(p: BigInt) -> (BigInt, BigInt, BigInt) {
            let fp = prime_field!(p);
            let coeff_a = generate_random_int_in_range(&BigInt::zero(), &(fp.p() - BigInt::one()));
            let coeff_a = FieldElement::new(&fp, Some(&coeff_a));

            let coeff_b = generate_random_int_in_range(&BigInt::zero(), &(fp.p() - BigInt::one()));
            let coeff_b = FieldElement::new(&fp, Some(&coeff_b));

            // d = 4a^3 + 27b^2
            let disc = FieldElement::new(&fp, Some(&BigInt::from(4u32)))
                * coeff_a.clone().pow(3u32)
                + FieldElement::new(&fp, Some(&BigInt::from(27u32))) * coeff_b.clone().pow(2u32);

            (coeff_a.0, coeff_b.0, disc.0)
        }

        loop {
            let (coeff_a, coeff_b, disc) = calc_curve_params(field.p());
            if disc != BigInt::zero() {
                let coeff_a = FieldElement::new(field, Some(&coeff_a));
                let coeff_b = FieldElement::new(field, Some(&coeff_b));
                let disc = FieldElement::new(field, Some(&disc));
                return Self {
                    coeff_a,
                    coeff_b,
                    disc,
                };
            }
        }
    }
}

pub trait EllipticCurve<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> {
    // An elliptic curve is identified by it's initial points.
    // new() returns the initial points
    fn new(
        x: Option<FieldElement<'a, P>>,
        y: Option<FieldElement<'a, P>>,
        c: CurveParams<'a, P>,
    ) -> Self;
    fn infinity_point(c: CurveParams<'a, P>) -> CurvePoint<'a, P>;
    fn slope(&self, other: &CurvePoint<'a, P>) -> Option<FieldElement<'a, P>>;
    fn underlying_field(&self) -> &P;
    fn get_coordinate_x(&self) -> Option<FieldElement<'a, P>>;
    fn get_coordinate_y(&self) -> Option<FieldElement<'a, P>>;
    fn get_curve_coeff_a(&self) -> Option<FieldElement<'a, P>>;
    fn get_curve_coeff_b(&self) -> Option<FieldElement<'a, P>>;
    fn get_curve_params(&self) -> Option<CurveParams<'a, P>>;
}

#[derive(Clone, Debug)]
pub struct CurvePoint<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> {
    x: Option<FieldElement<'a, P>>,
    y: Option<FieldElement<'a, P>>,
    curve_params: CurveParams<'a, P>,
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Display for CurvePoint<'a, P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "({}, \n    {})",
            self.get_coordinate_x().unwrap().0,
            self.get_coordinate_y().unwrap().0,
        )
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> PartialEq for CurvePoint<'a, P> {
    fn eq(&self, other: &Self) -> bool {
        match self.curve_params.coeff_a == other.curve_params.coeff_a
            && self.curve_params.coeff_b == other.curve_params.coeff_b
        {
            true => self.x == other.x && self.y == other.y,
            false => false,
        }
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> EcDisplay for CurvePoint<'a, P> {
    fn format_curve(&self) -> anyhow::Result<String> {
        let x = self.get_coordinate_x().clone().unwrap().0;
        let y = self.get_coordinate_y().clone().unwrap().0;
        let p = self.get_coordinate_x().unwrap().1.p();
        let a = self.get_curve_coeff_a().unwrap().0;
        let b = self.get_curve_coeff_b().unwrap().0;
        let y_square = format!("Y{}", Superscript(2));
        let x_cube = format!("X{}", Superscript(3));
        let res = format!(
            "  {} = {} + {}X \n\t + {} \n\t mod {}\n  Generator Point: \n    ({}, \n    {})\n",
            y_square, x_cube, a, b, p, x, y
        );
        Ok(res)
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> EllipticCurve<'a, P> for CurvePoint<'a, P> {
    fn new(
        x: Option<FieldElement<'a, P>>,
        y: Option<FieldElement<'a, P>>,
        c: CurveParams<'a, P>,
    ) -> Self {
        if x.is_some() && y.is_some() {
            Self {
                x,
                y,
                curve_params: c,
            }
        } else {
            loop {
                let x = generate_random_int_in_range(
                    &BigInt::zero(),
                    &(c.coeff_a.1.p() - BigInt::one()),
                );
                let x = FieldElement::new(c.coeff_a.1, Some(&x));
                let z = x.clone().pow(3u32) + c.coeff_a.clone() * x.clone() + c.coeff_b.clone();
                let q = modular_pow(
                    &z.0,
                    &((c.coeff_a.1.p() - BigInt::one()) / BigInt::from(2u32)),
                    &c.coeff_a.1.p(),
                );
                if q == BigInt::one() {
                    let exp = (c.coeff_a.1.p() + BigInt::one()) / 4;
                    let y = modular_pow(&z.0, &exp, &c.coeff_a.1.p());
                    let y = FieldElement::new(c.coeff_a.1, Some(&y));
                    return Self {
                        x: Some(x),
                        y: Some(y),
                        curve_params: c,
                    };
                }
            }
        }
    }

    fn infinity_point(c: CurveParams<'a, P>) -> Self {
        Self {
            x: None,
            y: None,
            curve_params: c,
        }
    }

    fn slope(&self, other: &CurvePoint<'a, P>) -> Option<FieldElement<'a, P>> {
        let coeff_a = self.curve_params.coeff_a.clone();
        let _coeff_b = self.curve_params.coeff_b.clone();

        let x1 = self.x.clone().unwrap();
        let y1 = self.y.clone().unwrap();

        if self.eq(other) {
            let two_y1 = FieldElement::new(x1.1, Some(&BigInt::from(2u32))) * y1;
            // (3x1^2 + a)*(2y1)^-1
            if let Some(inv) = two_y1.inv() {
                let slope = (FieldElement::new(x1.1, Some(&BigInt::from(3u32))) * x1.pow(2u32)
                    + coeff_a)
                    * inv;
                Some(slope)
            } else {
                None
            }
        } else {
            let x2 = other.x.clone().unwrap();
            let y2 = other.y.clone().unwrap();
            let x2_minus_x1 = x2 - x1;
            let y2_minus_y1 = y2 - y1;
            if let Some(inv) = x2_minus_x1.inv() {
                let slope = y2_minus_y1 * inv;
                Some(slope)
            } else {
                None
            }
        }
    }

    fn get_coordinate_x(&self) -> Option<FieldElement<'a, P>> {
        self.x.clone()
    }

    fn get_coordinate_y(&self) -> Option<FieldElement<'a, P>> {
        self.y.clone()
    }

    fn get_curve_coeff_a(&self) -> Option<FieldElement<'a, P>> {
        Some(self.curve_params.coeff_a.clone())
    }

    fn get_curve_coeff_b(&self) -> Option<FieldElement<'a, P>> {
        Some(self.curve_params.coeff_b.clone())
    }

    fn underlying_field(&self) -> &P {
        self.x.clone().unwrap().1
    }

    fn get_curve_params(&self) -> Option<CurveParams<'a, P>> {
        Some(self.curve_params.to_owned())
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Add<CurvePoint<'a, P>> for CurvePoint<'a, P> {
    type Output = CurvePoint<'a, P>;
    fn add(self, other: CurvePoint<'a, P>) -> CurvePoint<'a, P> {
        let infinity_point = Self::infinity_point(self.curve_params.clone());
        if !self.eq(&infinity_point) && other.eq(&infinity_point) {
            return self;
        } else if self.eq(&infinity_point) && !other.eq(&infinity_point) {
            return other;
        } else if self.eq(&infinity_point) && other.eq(&infinity_point) {
            return infinity_point;
        }

        if let Some(slope) = self.slope(&other) {
            let x3 = slope.clone().pow(2u32) - self.x.clone().unwrap() - other.x.unwrap();
            let y3 = slope * (self.x.unwrap() - x3.clone()) - self.y.unwrap();

            CurvePoint {
                x: Some(x3),
                y: Some(y3),
                curve_params: self.curve_params,
            }
        } else {
            Self::infinity_point(self.curve_params)
        }
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Neg for CurvePoint<'a, P> {
    type Output = CurvePoint<'a, P>;

    fn neg(self) -> Self::Output {
        CurvePoint {
            x: Some(self.x.unwrap()),
            y: Some(-self.y.unwrap()),
            curve_params: self.curve_params,
        }
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Mul<BigInt> for CurvePoint<'a, P> {
    type Output = CurvePoint<'a, P>;

    fn mul(self, rhs: BigInt) -> Self::Output {
        let mut q = CurvePoint::infinity_point(self.curve_params.clone());
        let (zero, one, _two) = (BigInt::from(0u64), BigInt::from(1u64), BigInt::from(2u64));
        let mut x = rhs.clone();
        let mut p = self.clone();
        while x > zero {
            // Extract the LSB from the exp.
            if &x & &one == one.clone() {
                q = q + p.clone();
            }

            p = p.clone() + p;

            // Division by 2 to get the next digit
            x /= 2;
        }

        q
    }
}

impl<'a, P: PrimeFieldConfig<BigInt> + Clone + Debug> Mul<CurvePoint<'a, P>> for BigInt {
    type Output = CurvePoint<'a, P>;

    fn mul(self, rhs: CurvePoint<'a, P>) -> Self::Output {
        rhs * self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ec_add_points_case1() {
        let q = prime_field!(BigInt::from(37u32));
        let coeff_a = FieldElement::new(&q, Some(&BigInt::from(-5i32)));
        let coeff_b = FieldElement::new(&q, Some(&BigInt::from(8i32)));
        // d = 4a^3 + 27b^2
        let disc = FieldElement::new(&q, Some(&BigInt::from(4u32))) * coeff_a.clone().pow(3u32)
            + FieldElement::new(&q, Some(&BigInt::from(27u32))) * coeff_b.clone().pow(2u32);
        let curve1_params = CurveParams {
            coeff_a,
            coeff_b,
            disc,
        };
        let _curve1 = CurvePoint::new(None, None, curve1_params.clone());
        let p1 = CurvePoint {
            x: Some(FieldElement::new(&q, Some(&BigInt::from(6u32)))),
            y: Some(FieldElement::new(&q, Some(&BigInt::from(3u32)))),
            curve_params: curve1_params.clone(),
        };
        let p2 = CurvePoint {
            x: Some(FieldElement::new(&q, Some(&BigInt::from(9u32)))),
            y: Some(FieldElement::new(&q, Some(&BigInt::from(10u32)))),
            curve_params: curve1_params.clone(),
        };
        let p3 = CurvePoint {
            x: Some(FieldElement::new(&q, Some(&BigInt::from(11u32)))),
            y: Some(FieldElement::new(&q, Some(&BigInt::from(10u32)))),
            curve_params: curve1_params.clone(),
        };
        let result = p1 + p2;
        assert_eq!(result, p3);
    }

    #[test]
    fn test_ec_add_points_case_2() {
        let q = prime_field!(BigInt::from(37u32));
        let coeff_a = FieldElement::new(&q, Some(&BigInt::from(-5i32)));
        let coeff_b = FieldElement::new(&q, Some(&BigInt::from(8i32)));
        // d = 4a^3 + 27b^2
        let disc = FieldElement::new(&q, Some(&BigInt::from(4u32))) * coeff_a.clone().pow(3u32)
            + FieldElement::new(&q, Some(&BigInt::from(27u32))) * coeff_b.clone().pow(2u32);
        let curve1_params = CurveParams {
            coeff_a,
            coeff_b,
            disc,
        };
        let _curve1 = CurvePoint::new(None, None, curve1_params.clone());
        let p1 = CurvePoint::infinity_point(curve1_params.clone());
        let p2 = CurvePoint {
            x: Some(FieldElement::new(&q, Some(&BigInt::from(9u32)))),
            y: Some(FieldElement::new(&q, Some(&BigInt::from(10u32)))),
            curve_params: curve1_params.clone(),
        };
        let result = p1 + p2.clone();
        assert_eq!(result, p2);
    }

    #[test]
    fn test_ec_add_points_case_3() {
        let q = prime_field!(BigInt::from(37u32));
        let coeff_a = FieldElement::new(&q, Some(&BigInt::from(-5i32)));
        let coeff_b = FieldElement::new(&q, Some(&BigInt::from(8i32)));
        // d = 4a^3 + 27b^2
        let disc = FieldElement::new(&q, Some(&BigInt::from(4u32))) * coeff_a.clone().pow(3u32)
            + FieldElement::new(&q, Some(&BigInt::from(27u32))) * coeff_b.clone().pow(2u32);
        let curve1_params = CurveParams {
            coeff_a,
            coeff_b,
            disc,
        };
        let _curve1 = CurvePoint::new(None, None, curve1_params.clone());
        let p2 = CurvePoint::infinity_point(curve1_params.clone());
        let p1 = CurvePoint {
            x: Some(FieldElement::new(&q, Some(&BigInt::from(9u32)))),
            y: Some(FieldElement::new(&q, Some(&BigInt::from(10u32)))),
            curve_params: curve1_params.clone(),
        };
        let result = p1.clone() + p2;
        assert_eq!(result, p1);
    }

    #[test]
    fn test_ec_add_points_case_none_none() {
        let q = prime_field!(BigInt::from(37u32));
        let coeff_a = FieldElement::new(&q, Some(&BigInt::from(-5i32)));
        let coeff_b = FieldElement::new(&q, Some(&BigInt::from(8i32)));
        // d = 4a^3 + 27b^2
        let disc = FieldElement::new(&q, Some(&BigInt::from(4u32))) * coeff_a.clone().pow(3u32)
            + FieldElement::new(&q, Some(&BigInt::from(27u32))) * coeff_b.clone().pow(2u32);
        let curve1_params = CurveParams {
            coeff_a,
            coeff_b,
            disc,
        };
        let _curve1 = CurvePoint::new(None, None, curve1_params.clone());
        let p1 = CurvePoint::infinity_point(curve1_params.clone());
        let p2 = CurvePoint::infinity_point(curve1_params.clone());
        let p3 = CurvePoint::infinity_point(curve1_params.clone());
        let result = p1 + p2;
        assert_eq!(result, p3);
    }

    #[test]
    fn test_ec_point_scalar_multiply_case1() {
        let q = prime_field!(BigInt::from(37u32));
        let coeff_a = FieldElement::new(&q, Some(&BigInt::from(-5i32)));
        let coeff_b = FieldElement::new(&q, Some(&BigInt::from(8i32)));
        // d = 4a^3 + 27b^2
        let disc = FieldElement::new(&q, Some(&BigInt::from(4u32))) * coeff_a.clone().pow(3u32)
            + FieldElement::new(&q, Some(&BigInt::from(27u32))) * coeff_b.clone().pow(2u32);
        let curve1_params = CurveParams {
            coeff_a,
            coeff_b,
            disc,
        };
        let _curve1 = CurvePoint::new(None, None, curve1_params.clone());
        let p1 = CurvePoint {
            x: Some(FieldElement::new(&q, Some(&BigInt::from(6u32)))),
            y: Some(FieldElement::new(&q, Some(&BigInt::from(3u32)))),
            curve_params: curve1_params.clone(),
        };
        let _p2 = CurvePoint {
            x: Some(FieldElement::new(&q, Some(&BigInt::from(9u32)))),
            y: Some(FieldElement::new(&q, Some(&BigInt::from(10u32)))),
            curve_params: curve1_params.clone(),
        };
        let p3 = CurvePoint {
            x: Some(FieldElement::new(&q, Some(&BigInt::from(16u32)))),
            y: Some(FieldElement::new(&q, Some(&BigInt::from(19u32)))),
            curve_params: curve1_params.clone(),
        };
        let result = p1.clone() * BigInt::from(5u32);
        assert_eq!(result, p3);

        let result = BigInt::from(5u32) * p1;
        assert_eq!(result, p3);
    }
}
