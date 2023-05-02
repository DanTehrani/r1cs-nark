pub mod acc_prover;
pub mod acc_verifier;

use crate::{utils::hadamard_prod, CurveAffineExt, MultiCommitGens, PRNG};
use ff::{Field, PrimeField, PrimeFieldBits};
use halo2curves::FieldExt;
use poseidon_transcript::transcript::PoseidonTranscript;
use zeroize::DefaultIsZeroes;

#[derive(Debug)]
pub struct HadamardInstance<C>(pub C, pub C, pub C)
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes;

#[derive(Debug, Clone)]
pub struct HadamardWitness<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub a_vec: Vec<C::ScalarExt>,
    pub b_vec: Vec<C::ScalarExt>,
    pub w1: C::ScalarExt,
    pub w2: C::ScalarExt,
    pub w3: C::ScalarExt,
}

pub struct HadamardAccumulator<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub qx: HadamardInstance<C>,
    pub qw: HadamardWitness<C>,
}

#[derive(Debug)]
pub struct HadamardAccProof<C>(Vec<C>)
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes;

#[cfg(test)]
mod tests {
    use super::{acc_prover::HadamardAccProver, acc_verifier::HadamardAccVerifier, *};
    use crate::predicates::hadamard::hadamard::HadamardProver;
    use crate::MultiCommitGens;
    use halo2curves::secq256k1::Secq256k1Affine;
    use halo2curves::CurveAffine;
    use poseidon_transcript::sponge::SpongeCurve;

    #[test]
    fn test_hadamard_accumulation() {
        type C = Secq256k1Affine;

        let n = 3; // Number of proofs to accumulate
        let l = 10; // Size of the vector

        let gens = MultiCommitGens::<C>::new(l, b"test_hadamard_accumulation");

        let mut a = vec![Vec::with_capacity(l); n];
        let mut b = vec![Vec::with_capacity(l); n];

        for i in 0..n {
            for j in 0..l {
                a[i].push(<C as CurveAffine>::ScalarExt::from_u128(i as u128));
                b[i].push(<C as CurveAffine>::ScalarExt::from_u128((n - i) as u128));
            }
        }

        let mut hadamard_prover = HadamardProver::new(gens.clone());
        let mut hadamard_instances = Vec::with_capacity(n);
        let mut hadamard_witnesses = Vec::with_capacity(n);

        for i in 0..n {
            let proof = hadamard_prover.prove(&a[i], &b[i]);
            let instance = HadamardInstance(proof.c1, proof.c2, proof.c3);
            let witness = HadamardWitness {
                a_vec: a[i].clone(),
                b_vec: b[i].clone(),
                w1: proof.w1,
                w2: proof.w2,
                w3: proof.w3,
            };

            hadamard_instances.push(instance);
            hadamard_witnesses.push(witness);
        }

        let acc_prover_transcript =
            PoseidonTranscript::new(b"test_hadamard_accumulation", SpongeCurve::K256);
        let mut acc_prover = HadamardAccProver::new(gens.clone(), acc_prover_transcript);

        let (acc, acc_proof) = acc_prover.prove_acc(&hadamard_instances, &hadamard_witnesses);

        let acc_verifier_transcript =
            PoseidonTranscript::new(b"test_hadamard_accumulation", SpongeCurve::K256);

        let mut acc_verifier = HadamardAccVerifier::new(acc_verifier_transcript);
        acc_verifier.verify(&acc.qx, &hadamard_instances, &acc_proof);
    }
}
