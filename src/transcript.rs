use super::r1cs::R1CS;
use crate::commitment::MultiCommitGens;
use crate::r1cs::{Pi1, Pi2, Proof};
use crate::utils::hadamard_prod;
use crate::CurveAffineExt;
use crate::PRNG;
use ff::{Field, PrimeField, PrimeFieldBits};
use poseidon_transcript::transcript::PoseidonTranscript;
use zeroize::DefaultIsZeroes;

pub struct Transcript<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes, {}

impl<C> Transcript<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    fn append_pi1() {}
}
