use super::utils::pi_1_to_hadamard_instance;
use crate::accumulation::hadamard::acc_prover::HadamardAccProver;
use crate::accumulation::hadamard::{HadamardAccProof, HadamardInstance, HadamardWitness};
use crate::r1cs::{Pi1, R1CSNARKProof, R1CS};
use crate::{CurveAffineExt, MultiCommitGens};
use ff::{Field, PrimeField, PrimeFieldBits};
use halo2curves::FieldExt;
use poseidon_transcript::sponge::SpongeCurve;
use poseidon_transcript::transcript::PoseidonTranscript;
use zeroize::DefaultIsZeroes;

#[derive(Debug)]
pub struct R1CSAccInstance<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub C_x: C,
    pub C_A: C,
    pub C_B: C,
    pub C_C: C,
    pub acc_HP_x: HadamardInstance<C>,
}

#[derive(Debug)]
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
    r1cs: R1CS<C>,
    gens: MultiCommitGens<C>,
    acc_transcript: PoseidonTranscript<C>,
    r1cs_transcript: PoseidonTranscript<C>,
}

#[derive(Debug)]
pub struct R1CSAccumulator<C>(R1CSAccInstance<C>, R1CSAccWitness<C>)
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
    C::Base: PrimeField<Repr = [u8; 32]>;

impl<C> R1CSAccProver<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
    C::Base: PrimeField<Repr = [u8; 32]>,
{
    pub fn new(
        r1cs: R1CS<C>,
        gens: MultiCommitGens<C>,
        acc_transcript: PoseidonTranscript<C>,
        r1cs_transcript: PoseidonTranscript<C>,
    ) -> Self {
        Self {
            r1cs,
            gens,
            r1cs_transcript,
            acc_transcript,
        }
    }

    pub fn prove_acc(
        &mut self,
        nizk_proofs: &[R1CSNARKProof<C>],
    ) -> (R1CSAccumulator<C>, HadamardAccProof<C>) {
        let n = nizk_proofs.len();
        let num_vars = self.r1cs.num_vars;
        let num_input = self.r1cs.num_input;

        let mut hadamard_instances = Vec::with_capacity(n);
        let mut hadamard_witnesses = Vec::with_capacity(n);
        let mut gammas = Vec::with_capacity(n);

        for i in 0..n {
            let r1cs_nizk_proof = &nizk_proofs[i];
            let pi_1 = &r1cs_nizk_proof.pi_1;
            let pi_2 = &r1cs_nizk_proof.pi_2;
            let public_input = &r1cs_nizk_proof.public_input;

            // TODO: Is this correct?
            self.r1cs_transcript
                .append_points(&[pi_1.C_A, pi_1.C_B, pi_1.C_C]);

            let gamma: C::ScalarExt = self.r1cs_transcript.squeeze(1)[0];
            gammas.push(gamma.clone());

            hadamard_instances.push(pi_1_to_hadamard_instance(&gamma, pi_1));
            self.r1cs_transcript.reset();

            let mut s_with_pub_inputs = Vec::with_capacity(num_vars + num_input);
            s_with_pub_inputs.extend_from_slice(&public_input);
            s_with_pub_inputs.extend_from_slice(&pi_2.s);

            let a_vec = self.r1cs.A.mul_vector(num_vars, &s_with_pub_inputs);
            let b_vec = self.r1cs.B.mul_vector(num_vars, &s_with_pub_inputs);
            let w1 = pi_2.sigma_A;
            let w2 = pi_2.sigma_B;
            let w3 = pi_2.sigma_O;

            hadamard_witnesses.push(HadamardWitness {
                a_vec,
                b_vec,
                w1,
                w2,
                w3,
            });
        }

        let hadamard_acc_prover_transcript =
            PoseidonTranscript::new(b"hadamard_prover", SpongeCurve::K256);
        let mut hadamard_acc_prover =
            HadamardAccProver::new(self.gens.clone(), hadamard_acc_prover_transcript);

        let (hadamard_acc, hadamard_acc_proof) =
            hadamard_acc_prover.prove_acc(&hadamard_instances, &hadamard_witnesses);

        // Step 5

        for hadamard_inst in hadamard_instances {
            self.acc_transcript.append_point(&hadamard_inst.0);
            self.acc_transcript.append_point(&hadamard_inst.1);
            self.acc_transcript.append_point(&hadamard_inst.2);
        }

        let beta: C::ScalarExt = self.acc_transcript.squeeze(1)[0];

        // Step 6

        let mut C_x: C::Curve = C::identity().into();
        let mut C_A: C::Curve = C::identity().into();
        let mut C_B: C::Curve = C::identity().into();
        let mut C_C: C::Curve = C::identity().into();

        let beta_pows = (0..n)
            .map(|i| beta.pow(&[i as u64, 0, 0, 0]))
            .collect::<Vec<C::ScalarExt>>();

        for (i, proof) in nizk_proofs.iter().enumerate() {
            C_x += self.gens.commit(&proof.public_input, &C::ScalarExt::zero()) * beta_pows[i];
            C_A += proof.pi_1.C_A + (proof.pi_1.C_A_prime * gammas[i]).into();
            C_B += proof.pi_1.C_B + (proof.pi_1.C_B_prime * gammas[i]).into();
            C_C += proof.pi_1.C_C + (proof.pi_1.C_C_prime * gammas[i]).into();
        }

        // Step 7

        let mut x = vec![C::ScalarExt::zero(); num_input];
        let mut s = vec![C::ScalarExt::zero(); num_vars];
        let mut sigma_A = C::ScalarExt::zero();
        let mut sigma_B = C::ScalarExt::zero();
        let mut sigma_C = C::ScalarExt::zero();

        for (i, proof) in nizk_proofs.iter().enumerate() {
            for k in 0..num_input {
                x[k] += proof.public_input[k] * beta_pows[i];
            }
            for k in 0..num_vars {
                s[k] += proof.pi_2.s[k] * beta_pows[i];
            }
            sigma_A += proof.pi_2.sigma_A * beta_pows[i];
            sigma_B += proof.pi_2.sigma_B * beta_pows[i];
            sigma_C += proof.pi_2.sigma_C * beta_pows[i];
        }

        let acc_witness = R1CSAccWitness {
            x,
            s,
            sigma_A,
            sigma_B,
            sigma_C,
            acc_HP_w: hadamard_acc.qw,
        };

        let acc_instance = R1CSAccInstance::<C> {
            C_x: C_x.into(),
            C_A: C_A.into(),
            C_B: C_B.into(),
            C_C: C_C.into(),
            acc_HP_x: hadamard_acc.qx,
        };

        let acc = R1CSAccumulator::<C>(acc_instance, acc_witness);

        (acc, hadamard_acc_proof)
    }
}
