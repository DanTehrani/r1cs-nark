use crate::CurveAffineExt;
use ff::{PrimeField, PrimeFieldBits};
use poseidon_transcript::sponge::{PoseidonSponge, SpongeCurve};
use rand_core::{OsRng, RngCore};
use zeroize::DefaultIsZeroes;

pub struct PRNG<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    sponge: PoseidonSponge<C::ScalarExt>,
}

impl<C> PRNG<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn new() -> Self {
        let mut sponge = PoseidonSponge::construct(b"r1cs-nark-prng", SpongeCurve::K256, None);
        let mut rng = OsRng;
        let mut bytes = [0u8; 32];
        rng.fill_bytes(&mut bytes);
        let random_scalar = C::ScalarExt::from_repr_vartime(bytes).unwrap();
        sponge.absorb(&[random_scalar]);

        Self { sponge }
    }

    pub fn squeeze(&mut self, length: usize) -> Vec<C::ScalarExt> {
        self.sponge.squeeze(length)
    }
}
