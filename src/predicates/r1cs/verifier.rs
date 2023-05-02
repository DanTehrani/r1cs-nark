use super::R1CS;
use crate::CurveAffineExt;
use crate::{r1cs::R1CSNARKProof, utils::hadamard_prod, MultiCommitGens};
use ff::{PrimeField, PrimeFieldBits};
use halo2curves::group::Curve;
use halo2curves::FieldExt;
pub use poseidon_transcript::transcript::PoseidonTranscript;
use zeroize::DefaultIsZeroes;

pub struct Verifier<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
    C::Base: FieldExt<Repr = [u8; 32]>,
{
    transcript: PoseidonTranscript<C>,
    pub r1cs: R1CS<C>,
    pub comm_gens: MultiCommitGens<C>,
}

impl<C> Verifier<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
    C::Base: FieldExt<Repr = [u8; 32]>,
{
    pub fn new(
        r1cs: R1CS<C>,
        transcript: PoseidonTranscript<C>,
        comm_gens: MultiCommitGens<C>,
    ) -> Self {
        Self {
            r1cs,
            transcript,
            comm_gens,
        }
    }

    pub fn verify(&mut self, proof: &R1CSNARKProof<C>, public_input: &Vec<C::ScalarExt>) {
        // Verify following the steps described in Section 8.1 of

        let num_cons = self.r1cs.num_cons;
        let num_vars = self.r1cs.num_vars;

        // Step 1

        let pi_1 = &proof.pi_1;
        let pi_2 = &proof.pi_2;

        // Step 2

        self.transcript.append_points(&[
            pi_1.C_A,
            pi_1.C_B,
            pi_1.C_C,
            pi_1.C_A_prime,
            pi_1.C_B_prime,
            pi_1.C_C_prime,
            pi_1.C_1,
            pi_1.C_2,
        ]);

        let gamma: C::ScalarExt = self.transcript.squeeze(1)[0];

        // Step 3

        let mut s_with_pub_input = Vec::with_capacity(num_vars + public_input.len());
        s_with_pub_input.extend_from_slice(&public_input);
        s_with_pub_input.extend_from_slice(&pi_2.s);

        assert_eq!(s_with_pub_input.len(), num_vars + public_input.len());

        let s_A = self.r1cs.A.mul_vector(num_cons, &s_with_pub_input);
        let s_B = self.r1cs.B.mul_vector(num_cons, &s_with_pub_input);
        let s_C = self.r1cs.C.mul_vector(num_cons, &s_with_pub_input);

        // Step 4

        let comm_s_A: C::Curve = self.comm_gens.commit(&s_A, &pi_2.sigma_A).into();
        let comm_s_B: C::Curve = self.comm_gens.commit(&s_B, &pi_2.sigma_B).into();
        let comm_s_C: C::Curve = self.comm_gens.commit(&s_C, &pi_2.sigma_C).into();

        assert_eq!(comm_s_A, pi_1.C_A + (pi_1.C_A_prime * gamma).into());
        assert_eq!(comm_s_B, pi_1.C_B + (pi_1.C_B_prime * gamma).into());
        assert_eq!(comm_s_C, pi_1.C_C + (pi_1.C_C_prime * gamma).into());

        // Step 5

        let comm_s_A_s_B: C::Curve = self
            .comm_gens
            .commit(&hadamard_prod::<C>(&s_A, &s_B), &pi_2.sigma_O)
            .into();

        println!("comm_s_A_s_B: {:?}", comm_s_A_s_B.to_affine());

        let C_C: C::Curve = pi_1.C_C.into();
        let C_1_gamma: C::Curve = (pi_1.C_1 * gamma).into();
        let C_2_gamma_squared: C::Curve = (pi_1.C_2 * gamma * gamma).into();

        let lhs = C_C + C_1_gamma + C_2_gamma_squared;

        println!("lhs: {:?}", lhs.to_affine());

        assert_eq!(comm_s_A_s_B, C_C + C_1_gamma + C_2_gamma_squared);
    }
}
