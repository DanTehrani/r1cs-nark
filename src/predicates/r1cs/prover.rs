use super::r1cs::R1CS;
use crate::commitment::MultiCommitGens;
use crate::r1cs::{Pi1, Pi2, Proof};
use crate::utils::hadamard_prod;
use crate::CurveAffineExt;
use crate::PRNG;
use ff::{Field, PrimeField, PrimeFieldBits};
use poseidon_transcript::transcript::PoseidonTranscript;
use zeroize::DefaultIsZeroes;

pub struct Prover<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub r1cs: R1CS<C>,
    pub transcript: PoseidonTranscript<C>,
    pub prng: PRNG<C>,
    pub comm_gens: MultiCommitGens<C>,
}

impl<C> Prover<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn new(
        r1cs: R1CS<C>,
        transcript: PoseidonTranscript<C>,
        comm_gens: MultiCommitGens<C>,
    ) -> Self {
        let prng = PRNG::new();
        Self {
            r1cs,
            transcript,
            prng,
            comm_gens,
        }
    }

    pub fn prove(
        &mut self,
        witness: &Vec<C::ScalarExt>,
        public_input: &Vec<C::ScalarExt>,
    ) -> Proof<C> {
        let num_cons = self.r1cs.num_cons;
        let num_vars = self.r1cs.num_vars;

        // Prove following the steps described in Section 8.1 of

        // Step 1

        let mut z = Vec::with_capacity(witness.len() + public_input.len());
        z.extend(public_input);
        z.extend(witness);

        // Step 2

        let r = self.prng.squeeze(num_vars);

        // Step 3

        let z_A = self.r1cs.A.mul_vector(num_cons, &z);
        let z_B = self.r1cs.B.mul_vector(num_cons, &z);
        let z_C = self.r1cs.C.mul_vector(num_cons, &z);

        let mut padded_r = vec![C::ScalarExt::zero(); public_input.len()];
        padded_r.extend_from_slice(&r);

        let r_A = self.r1cs.A.mul_vector(num_cons, &padded_r);
        let r_B = self.r1cs.B.mul_vector(num_cons, &padded_r);
        let r_C = self.r1cs.C.mul_vector(num_cons, &padded_r);

        // Step 4

        let w = self.prng.squeeze(8);

        let w_A = w[0];
        let w_B = w[1];
        let w_C = w[2];

        let C_A = self.comm_gens.commit(&z_A, &w_A);
        let C_B = self.comm_gens.commit(&z_B, &w_B);
        let C_C = self.comm_gens.commit(&z_C, &w_C);

        let w_A_prime = w[3];
        let w_B_prime = w[4];
        let w_C_prime = w[5];

        let C_A_prime = self.comm_gens.commit(&r_A, &w_A_prime);
        let C_B_prime = self.comm_gens.commit(&r_B, &w_B_prime);
        let C_C_prime = self.comm_gens.commit(&r_C, &w_C_prime);

        let w_1 = w[6];
        let w_2 = w[7];

        // Step 5

        let z_A_r_B = hadamard_prod::<C>(&z_A, &r_B);
        let z_B_r_A = hadamard_prod::<C>(&z_B, &r_A);

        let cross_term_1: Vec<C::ScalarExt> = z_A_r_B
            .iter()
            .zip(z_B_r_A.iter())
            .map(|(a, b)| *a + *b)
            .collect();

        let cross_term_2: Vec<C::ScalarExt> = hadamard_prod::<C>(&r_A, &r_B);

        let C_1 = self.comm_gens.commit(&cross_term_1, &w_1);
        let C_2 = self.comm_gens.commit(&cross_term_2, &w_2);

        self.transcript
            .append_points(&[C_A, C_B, C_C, C_A_prime, C_B_prime, C_C_prime, C_1, C_2]);

        // Step 6

        let pi_1 = Pi1 {
            C_A,
            C_B,
            C_C,
            C_A_prime,
            C_B_prime,
            C_C_prime,
            C_1,
            C_2,
        };

        // Step 7

        let gamma = self.transcript.squeeze(1)[0];

        // Step 8

        let s = witness
            .iter()
            .zip(r.iter())
            .map(|(w_i, r_i)| *w_i + *r_i * gamma)
            .collect::<Vec<C::ScalarExt>>();

        assert_eq!(s.len(), witness.len());

        // Step 9

        let sigma_A = w_A + gamma * w_A_prime;
        let sigma_B = w_B + gamma * w_B_prime;
        let sigma_C = w_C + gamma * w_C_prime;

        // Step 10

        let sigma_O = w_C + (gamma * w_1) + (gamma * gamma * w_2);

        // Step 11

        let pi_2 = Pi2 {
            sigma_A,
            sigma_B,
            sigma_C,
            sigma_O,
            s,
        };

        // Step 12

        Proof { pi_1, pi_2 }
    }
}
