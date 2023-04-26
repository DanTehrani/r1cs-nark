use crate::CurveAffineExt;
use digest::{ExtendableOutput, Input};
use ff::{PrimeField, PrimeFieldBits};
use multiexp::multiexp;
use sha3::Shake256;
use std::io::Read;
use zeroize::DefaultIsZeroes;

#[derive(Debug)]
pub struct MultiCommitGens<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub G: Vec<C>,
    pub h: C,
}

impl<C> MultiCommitGens<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn new(n: usize, label: &[u8]) -> Self {
        let mut shake = Shake256::default();
        shake.input(label);
        shake.input(C::generator().to_bytes());

        let mut reader = shake.xof_result();
        let mut gens: Vec<C> = Vec::new();
        let mut uniform_bytes = [0u8; 128];
        for _ in 0..n + 1 {
            reader.read_exact(&mut uniform_bytes).unwrap();
            // TODO: Curve point from random bytes
            // Unsafe!
            gens.push(C::generator());
        }

        MultiCommitGens {
            G: gens[..n].to_vec(),
            h: gens[n],
        }
    }

    pub fn clone(&self) -> Self {
        Self {
            h: self.h,
            G: self.G.clone(),
        }
    }

    pub fn commit(&self, a: &Vec<C::ScalarExt>, blinder: &C::ScalarExt) -> C {
        assert_eq!(self.G.len(), a.len());

        let mut com = self.h * blinder;

        let pairs = self
            .G
            .iter()
            .enumerate()
            .map(|(i, g)| (a[i], (*g).into()))
            .collect::<Vec<(C::ScalarExt, C::Curve)>>();

        let com = multiexp(&pairs);
        com.into()
    }
}
