use crate::accumulation::hadamard::{
    HadamardAccProof, HadamardAccumulator, HadamardInstance, HadamardWitness,
};
use crate::{CurveAffineExt, MultiCommitGens};
use ff::{Field, PrimeField, PrimeFieldBits};
use halo2curves::FieldExt;
use poseidon_transcript::sponge::SpongeCurve;
use poseidon_transcript::transcript::PoseidonTranscript;
use zeroize::DefaultIsZeroes;

pub struct HadamardAccVerifier<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    transcript: PoseidonTranscript<C>,
}

impl<C> HadamardAccVerifier<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn new(transcript: PoseidonTranscript<C>) -> Self {
        Self { transcript }
    }

    pub fn verify(
        &mut self,
        acc_instance: &HadamardInstance<C>,
        accumulated_instances: &[HadamardInstance<C>],
        proof: &HadamardAccProof<C>,
    ) {
        let n = accumulated_instances.len();

        // Absorb the accumulator instances
        for acc_inst in accumulated_instances {
            self.transcript.append_point(&acc_inst.0);
            self.transcript.append_point(&acc_inst.1);
            self.transcript.append_point(&acc_inst.2);
        }

        let mu: C::ScalarExt = self.transcript.squeeze(1)[0];
        let nu: C::ScalarExt = self.transcript.squeeze(1)[0];

        let mut mu_powers = vec![];
        for i in 0..n {
            mu_powers.push(mu.pow(&[i as u64, 0, 0, 0]));
        }

        let mut nu_powers = vec![];
        for i in 0..n {
            nu_powers.push(nu.pow(&[i as u64, 0, 0, 0]));
        }

        let mut expected_c1: C::Curve = C::identity().into();
        let mut expected_c2: C::Curve = C::identity().into();
        for i in 0..n {
            expected_c1 += accumulated_instances[i].0 * mu_powers[i] * nu_powers[i];
            expected_c2 += accumulated_instances[i].1 * nu_powers[nu_powers.len() - 1 - i];
        }

        let mut expected_c3_1: C::Curve = C::identity().into();
        let mut expected_c3_2: C::Curve = C::identity().into();
        let mut expected_c3_3: C::Curve = C::identity().into();

        for i in 0..(n - 1) {
            expected_c3_1 += proof.0[i] * nu_powers[i];
        }

        for i in 0..accumulated_instances.len() {
            expected_c3_2 += accumulated_instances[i].2 * mu_powers[i];
        }
        expected_c3_2 *= nu_powers[n - 1];

        for i in 0..(n - 1) {
            expected_c3_3 += proof.0[i + n - 1] * nu.pow(&[(n + i) as u64, 0, 0, 0]);
        }

        let expected_c3 = expected_c3_1 + expected_c3_2 + expected_c3_3;

        assert_eq!(expected_c1, acc_instance.0.into());
        assert_eq!(expected_c2, acc_instance.1.into());
        assert_eq!(expected_c3, acc_instance.2.into());
    }
}
