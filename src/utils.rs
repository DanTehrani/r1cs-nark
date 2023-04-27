use crate::CurveAffineExt;
use ff::{Field, PrimeField, PrimeFieldBits};
use zeroize::DefaultIsZeroes;

pub fn hadamard_prod<C>(a: &Vec<C::ScalarExt>, b: &Vec<C::ScalarExt>) -> Vec<C::ScalarExt>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    assert_eq!(a.len(), b.len());
    let mut result = vec![C::ScalarExt::zero(); a.len()];
    for i in 0..a.len() {
        result[i] = a[i] * b[i];
    }
    result
}
