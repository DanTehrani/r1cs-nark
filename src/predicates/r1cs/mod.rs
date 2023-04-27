#![allow(non_snake_case)]
mod prover;
mod r1cs;
mod verifier;

pub use prover::Prover;
pub use r1cs::R1CS;
pub use verifier::Verifier;

use crate::CurveAffineExt;
use ff::{PrimeField, PrimeFieldBits};
use zeroize::DefaultIsZeroes;

pub struct Pi1<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    C_A: C,
    C_B: C,
    C_C: C,
    C_A_prime: C,
    C_B_prime: C,
    C_C_prime: C,
    C_1: C,
    C_2: C,
}

pub struct Pi2<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    s: Vec<C::ScalarExt>,
    sigma_A: C::ScalarExt,
    sigma_B: C::ScalarExt,
    sigma_C: C::ScalarExt,
    sigma_O: C::ScalarExt,
}

pub struct Proof<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub pi_1: Pi1<C>,
    pub pi_2: Pi2<C>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commitment::MultiCommitGens;
    use crate::r1cs::{Prover, Verifier};
    use halo2curves::secq256k1::Secq256k1Affine;
    use poseidon_transcript::sponge::SpongeCurve;
    use poseidon_transcript::transcript::PoseidonTranscript;

    #[test]
    pub fn test_r1cs_nizk() {
        let num_cons = 10;
        let num_vars = 10;
        let num_input = 5;

        type C = Secq256k1Affine;

        let prover_transcript = PoseidonTranscript::new(b"test-r1cs-nizk", SpongeCurve::K256);
        let r1cs = R1CS::<C>::produce_synthetic_r1cs(num_cons, num_vars, num_input);

        let gens = MultiCommitGens::<C>::new(num_cons, b"r1cs-nark");
        let mut prover = Prover::new(r1cs.clone(), prover_transcript, gens.clone());

        let proof = prover.prove(&r1cs.witness, &r1cs.public_input);

        let verifier_transcript = PoseidonTranscript::new(b"test-r1cs-nizk", SpongeCurve::K256);
        let mut verifier = Verifier::new(r1cs.clone(), verifier_transcript, gens.clone());

        verifier.verify(&proof, &r1cs.public_input);
    }
}
