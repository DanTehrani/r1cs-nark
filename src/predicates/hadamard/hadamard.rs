use crate::{utils::hadamard_prod, CurveAffineExt, MultiCommitGens, PRNG};
use ff::{PrimeField, PrimeFieldBits};
use zeroize::DefaultIsZeroes;

pub struct HadamardProof<C: CurveAffineExt> {
    pub c1: C,
    pub c2: C,
    pub c3: C,
    pub w1: C::ScalarExt,
    pub w2: C::ScalarExt,
    pub w3: C::ScalarExt,
}

pub struct HadamardProver<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    gens: MultiCommitGens<C>,
    prng: PRNG<C>,
}

impl<C> HadamardProver<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn new(gens: MultiCommitGens<C>) -> Self {
        Self {
            gens,
            prng: PRNG::new(),
        }
    }

    pub fn prove(&mut self, a: &[C::ScalarExt], b: &[C::ScalarExt]) -> HadamardProof<C> {
        let c = hadamard_prod::<C>(&a, &b);

        let w = self.prng.squeeze(3);
        let w1 = w[0];
        let w2 = w[1];
        let w3 = w[2];

        let c1 = self.gens.commit(&a, &w1);
        let c2 = self.gens.commit(&b, &w2);
        let c3 = self.gens.commit(&c, &w3);

        HadamardProof {
            c1,
            c2,
            c3,
            w1: w1,
            w2: w2,
            w3: w3,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use halo2curves::{secq256k1::Secq256k1Affine, CurveAffine, FieldExt};

    // We don't use a HadamardProof verifier in any of the protocol we implement,
    // so we just write a simple test here to make sure the proof verifies.
    #[test]
    fn test_hadamard_prove_and_verify() {
        type C = Secq256k1Affine;

        let n = 10;
        let gens: MultiCommitGens<_> = MultiCommitGens::new(n, b"test_hadamard_proof_verify");
        let mut hadamard_prover = HadamardProver::<C>::new(gens.clone());

        let mut a = Vec::with_capacity(n);
        let mut b = Vec::with_capacity(n);

        for i in 0..n {
            a.push(<C as CurveAffine>::ScalarExt::from_u128(i as u128));
            b.push(<C as CurveAffine>::ScalarExt::from_u128((n - i) as u128));
        }

        let proof = hadamard_prover.prove(&a, &b);

        assert_eq!(gens.commit(&a, &proof.w1), proof.c1);
        assert_eq!(gens.commit(&b, &proof.w2), proof.c2);
        assert_eq!(
            gens.commit(&hadamard_prod::<C>(&a, &b), &proof.w3),
            proof.c3
        );
    }
}
