//! Naive implementation of the KZG Scheme over the BLS12_381 curve, zero optimizations
//! The perdersen version of the scheme is not implemented

use anyhow::anyhow;
use ark_bls12_381::Fr;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, Group};
use ark_ff::{BigInteger128, BigInteger256, Field, One, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{
    DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial as PolynomialT,
};
use std::ops::{Div, Sub};

pub type G1Affine = ark_bls12_381::g1::G1Affine;
pub type G1Projective = ark_bls12_381::g1::G1Projective;
pub type G2Affine = ark_bls12_381::g2::G2Affine;
pub type G2Projective = ark_bls12_381::g2::G2Projective;
pub type Polynomial = DensePolynomial<Fr>;

#[derive(Debug)]
pub struct Proof {
    pub commitment: Commitment,
    pub evaluation: Evaluation,
    pub x: Fr,
}

#[derive(Debug)]
pub struct BatchProof {
    pub commitment: Commitment,
    pub evaluations: Vec<Evaluation>,
}

pub type Evaluation = Fr;

pub type Commitment = G1Projective;

#[derive(Debug)]
pub struct KZG {
    pub powers_of_g1: Vec<G1Projective>,
    pub powers_of_g2: Vec<G2Projective>,
    pub commitment: Commitment,
}

impl KZG {
    pub fn setup(random_scalar: Fr, poly: &Polynomial) -> Self {
        let g1 = G1Affine::new(
            ark_bls12_381::g1::G1_GENERATOR_X,
            ark_bls12_381::g1::G1_GENERATOR_Y,
        );
        let g2 = G2Affine::new(
            ark_bls12_381::g2::G2_GENERATOR_X,
            ark_bls12_381::g2::G2_GENERATOR_Y,
        );
        let deg = poly.degree();
        let mut powers_of_g1 = vec![];
        let mut powers_of_g2 = vec![];
        for i in 0..=deg {
            let g1_i = g1.mul_bigint(BigInteger256::from(
                random_scalar.pow(BigInteger128::from(i as u64)),
            ));
            let g2_i = g2.mul_bigint(BigInteger256::from(
                random_scalar.pow(BigInteger128::from(i as u64)),
            ));
            powers_of_g1.push(g1_i);
            powers_of_g2.push(g2_i)
        }

        let commitment =
            Self::commit(&powers_of_g1, &poly).expect("Failed to commit to the given polynomial");

        Self {
            powers_of_g1,
            powers_of_g2,
            commitment,
        }
    }

    fn commit(powers_of_g: &[G1Projective], poly: &Polynomial) -> anyhow::Result<Commitment> {
        if powers_of_g.len() < poly.degree() {
            return Err(anyhow!(
                "Degree of polynomial is greater than powers of G supplied"
            ));
        }
        let mut commitment: G1Projective = G1Projective::zero();
        for (i, coeff) in poly.coeffs.iter().enumerate() {
            let g_i = powers_of_g[i];
            commitment = commitment + g_i.mul_bigint(BigInteger256::from(*coeff))
        }
        Ok(commitment)
    }

    pub fn point_evaluation_proof(&self, poly: &Polynomial, x: Fr) -> anyhow::Result<Proof> {
        let evaluation = poly.evaluate(&x);
        let dividend = poly.sub(&DensePolynomial::from_coefficients_vec(vec![evaluation]));
        let divisor = DensePolynomial::from_coefficients_vec(vec![-x, Fr::one()]);
        let quotient = dividend.div(&divisor);
        let commitment = Self::commit(&self.powers_of_g1, &quotient)?;
        Ok(Proof {
            commitment,
            evaluation,
            x,
        })
    }

    fn batch_evaluation_proof(&self, poly: Polynomial, values: Vec<Fr>) -> BatchProof {
        todo!()
    }

    /// e(CQ,a⋅G2−b⋅G2)=?e(CF−G1⋅c,G2)
    pub fn verify_single_eval(&self, proof: Proof) -> bool {
        let g1: G1Projective = G1Affine::new(
            ark_bls12_381::g1::G1_GENERATOR_X,
            ark_bls12_381::g1::G1_GENERATOR_Y,
        )
        .into();
        let g2: G2Projective = G2Affine::new(
            ark_bls12_381::g2::G2_GENERATOR_X,
            ark_bls12_381::g2::G2_GENERATOR_Y,
        )
        .into();

        // a.G2
        let a_g2 = self.powers_of_g2[1].clone();
        // b.G2
        let b_g2 = g2.mul_bigint(BigInteger256::from(proof.x));

        // G1.c
        let c_g1 = g1.mul_bigint(BigInteger256::from(proof.evaluation));

        let p_1 = ark_bls12_381::Bls12_381::pairing(proof.commitment, a_g2 - b_g2);
        let p_2 = ark_bls12_381::Bls12_381::pairing(self.commitment - c_g1, g2);
        p_1 == p_2
    }

    fn verify_batch_eval(&self, batch_proof: BatchProof) -> bool {
        todo!()
    }
}

pub fn lagrange_interpolation(deg: usize, evals: Vec<Evaluation>) -> anyhow::Result<Polynomial> {
    todo!()
}

pub fn fft(num_coeffs: usize, mut evals: Vec<Evaluation>) -> anyhow::Result<Polynomial> {
    let domain = GeneralEvaluationDomain::<Fr>::new(num_coeffs)
        .ok_or_else(|| anyhow!("Polynomial interpolation failed"))?;
    domain.ifft_in_place(&mut evals);
    Ok(Polynomial { coeffs: evals })
}

#[cfg(test)]
mod test {
    use crate::{Polynomial, KZG};
    use ark_bls12_381::Fr;
    use ark_poly::DenseUVPolynomial;
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn test_single_proof_eval() {
        let mut test_rng = test_rng();
        let poly = Polynomial::rand(100, &mut test_rng);

        let sk = Fr::rand(&mut test_rng);

        let kzg = KZG::setup(sk, &poly);
        let x = Fr::rand(&mut test_rng);
        let proof = kzg.point_evaluation_proof(&poly, x).unwrap();

        assert!(kzg.verify_single_eval(proof))
    }
}
