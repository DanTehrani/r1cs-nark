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
    C::ScalarExt: DefaultIsZeroes,
{
    fn new(gens: MultiCommitGens<C>) -> Self {
        Self {
            gens,
            prng: PRNG::new(),
        }
    }

    pub fn prove(&mut self, a: Vec<C::ScalarExt>, b: Vec<C::ScalarExt>) -> HadamardProof<C> {
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
