use crate::accumulation::hadamard::acc_prover::HadamardAccProver;
use crate::accumulation::hadamard::{HadamardAccProof, HadamardInstance, HadamardWitness};
use crate::r1cs::{Pi1, R1CSNARKProof, R1CS};
use crate::{CurveAffineExt, MultiCommitGens};
use ff::{Field, PrimeField, PrimeFieldBits};
use halo2curves::FieldExt;
use poseidon_transcript::sponge::SpongeCurve;
use poseidon_transcript::transcript::PoseidonTranscript;
use zeroize::DefaultIsZeroes;

pub fn pi_1_to_hadamard_instance<C>(gamma: &C::ScalarExt, pi_1: &Pi1<C>) -> HadamardInstance<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    let C_C: C::Curve = pi_1.C_C.into();
    let C_1_gamma: C::Curve = (pi_1.C_1 * gamma).into();
    let C_2_gamma_squared: C::Curve = (pi_1.C_2 * (*gamma * *gamma)).into();

    HadamardInstance::<C>(
        (pi_1.C_A + (pi_1.C_A_prime * gamma).into()).into(),
        (pi_1.C_B + (pi_1.C_B_prime * gamma).into()).into(),
        (C_C + C_1_gamma + C_2_gamma_squared).into(),
    )
}
