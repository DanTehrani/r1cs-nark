use criterion::{black_box, criterion_group, criterion_main, Criterion};
use halo2curves::secq256k1::Secq256k1Affine;
use poseidon_transcript::sponge::SpongeCurve;
use poseidon_transcript::transcript::PoseidonTranscript;
use r1cs_nark::r1cs::{R1CSNARKProver, R1CS};
use r1cs_nark::MultiCommitGens;

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("r1cs-prove");
    type C = Secq256k1Affine;

    for s in &[10, 13, 15] {
        let num_cons = 2usize.pow(*s);
        let num_vars = num_cons;
        let num_input = 10;

        let prover_transcript = PoseidonTranscript::<C>::new(b"test-prove", SpongeCurve::K256);
        let r1cs = R1CS::produce_synthetic_r1cs(num_cons, num_vars, num_input);
        let gens = MultiCommitGens::new(num_cons, b"r1cs-nark");

        let mut prover = R1CSNARKProver::new(r1cs.clone(), prover_transcript, gens);

        let witness = &r1cs.witness;
        let public_input = &r1cs.public_input;

        let name = format!("Prove {} constraints", num_cons);
        group.bench_function(name, move |b| {
            b.iter(|| {
                prover.prove(black_box(witness), black_box(public_input));
            });
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
