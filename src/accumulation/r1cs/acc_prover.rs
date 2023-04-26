use crate::accumulation::hadamard::acc_prover::{HadamardInstance, HadamardWitness};
use crate::{utils::hadamard_prod, CurveAffineExt, MultiCommitGens, PRNG};
use ff::{Field, PrimeField, PrimeFieldBits};
use halo2curves::FieldExt;
use poseidon_transcript::sponge::PoseidonSponge;
use zeroize::DefaultIsZeroes;

pub struct R1CSAccInstance<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub C_x: C,
    pub C_A: C,
    pub C_B: C,
    pub C_C: C,
    pub acc_HP_x: HadamardInstance<C>,
}

pub struct R1CSAccWitness<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub x: Vec<C::ScalarExt>,
    pub s: Vec<C::ScalarExt>,
    pub sigma_A: C::ScalarExt,
    pub sigma_B: C::ScalarExt,
    pub sigma_C: C::ScalarExt,
    pub acc_HP_w: HadamardWitness<C>,
}

pub struct R1CSAccProver<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    gens: MultiCommitGens<C>,
    sponge: PoseidonSponge<C::ScalarExt>,
}

impl<C> R1CSAccProver<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn new(gens: MultiCommitGens<C>, sponge: PoseidonSponge<C::ScalarExt>) -> Self {
        Self {
            gens,
            sponge,
            prng: PRNG::new(),
        }
    }

    pub fn prove_acc(
        &mut self,
        qx: Vec<HadamardInstance<C>>,
        qw: Vec<HadamardWitness<C>>,
    ) -> (R1CSAccumulator<C>, R1CSAccProof<C>) {
    }
}
