use crate::accumulation::hadamard::acc_verifier::HadamardAccVerifier;
use crate::accumulation::hadamard::{
    HadamardAccProof, HadamardAccumulator, HadamardInstance, HadamardWitness,
};
use crate::{CurveAffineExt, MultiCommitGens};
use ff::{Field, PrimeField, PrimeFieldBits};
use halo2curves::FieldExt;
use poseidon_transcript::sponge::SpongeCurve;
use poseidon_transcript::transcript::PoseidonTranscript;
use zeroize::DefaultIsZeroes;

use super::acc_prover::{R1CSAccInstance, R1CSAccumulator};
use super::utils::pi_1_to_hadamard_instance;

pub struct R1CSAccVerifier<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    acc_transcript: PoseidonTranscript<C>,
    r1cs_transcript: PoseidonTranscript<C>,
}

impl<C> R1CSAccVerifier<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn new(
        acc_transcript: PoseidonTranscript<C>,
        r1cs_transcript: PoseidonTranscript<C>,
    ) -> Self {
        Self {
            acc_transcript,
            r1cs_transcript,
        }
    }

    pub fn verify(
        &mut self,
        acc: &R1CSAccumulator<C>,
        accumulated_instances: &[R1CSAccInstance<C>],
        proof: &HadamardAccProof<C>,
    ) {
        let hadamard_acc_verifier = HadamardAccVerifier::new(self.acc_transcript.clone());
        let n = accumulated_instances.len();

        let mut gammas = Vec::with_capacity(n);

        //        let mut hadamard_instances = Vec::with_capacity(n);
        //        let mut hadamard_witnesses = Vec::with_capacity(n);

        for i in 0..n {
            let acc_instance = &accumulated_instances[i];

            // TODO: Is this correct?
            self.r1cs_transcript.append_points(&[
                acc_instance.C_A,
                acc_instance.C_B,
                acc_instance.C_C,
            ]);

            let gamma: C::ScalarExt = self.r1cs_transcript.squeeze(1)[0];
            gammas.push(gamma.clone());

            self.r1cs_transcript.reset();
        }
        // TBD
    }
}
