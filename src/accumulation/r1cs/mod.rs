mod acc_prover;
mod acc_verifier;
mod utils;

pub use acc_prover::R1CSAccProver;

#[cfg(test)]
mod tests {
    use crate::{
        r1cs::{R1CSNARKProver, R1CS},
        MultiCommitGens,
    };
    use halo2curves::secq256k1::Secq256k1Affine;
    use poseidon_transcript::{sponge::SpongeCurve, transcript::PoseidonTranscript};

    use super::*;

    #[test]
    fn test_r1cs_accumulation() {
        let num_cons = 10;
        let num_vars = 10;
        let num_input = 10;
        type C = Secq256k1Affine;

        let gens = MultiCommitGens::new(num_vars, b"test_r1cs_accumulation");

        let r1cs = R1CS::produce_synthetic_r1cs(num_cons, num_vars, num_input);

        let r1cs_nizk_prover_transcript =
            PoseidonTranscript::new(b"test_r1cs_accumulation", SpongeCurve::K256);
        let mut r1cs_nizk_prover =
            R1CSNARKProver::new(r1cs.clone(), r1cs_nizk_prover_transcript, gens.clone());

        // Generate the R1CS NARK proofs

        let n = 3; // Number of r1cs instances to accumulate

        let mut r1cs_nizk_proofs = Vec::with_capacity(n);
        for _ in 0..n {
            r1cs_nizk_proofs.push(r1cs_nizk_prover.prove(&r1cs.witness, &r1cs.public_input));
        }

        // Generate the accumulation proof

        let acc_prover_transcript =
            PoseidonTranscript::new(b"test_r1cs_accumulation", SpongeCurve::K256);

        let r1cs_nizk_prover_transcript =
            PoseidonTranscript::new(b"test_r1cs_accumulation", SpongeCurve::K256);

        let mut r1cs_acc_prover = R1CSAccProver::<C>::new(
            r1cs,
            gens.clone(),
            acc_prover_transcript,
            r1cs_nizk_prover_transcript,
        );

        let (accumulator, acc_proof) = r1cs_acc_prover.prove_acc(&r1cs_nizk_proofs);
        println!("Accumulator: {:?}", accumulator);
        println!("Accumulation proof: {:?}", acc_proof);
    }
}
