use crate::CurveAffineExt;
use ff::{Field, PrimeField, PrimeFieldBits};
use halo2curves::FieldExt;
use zeroize::DefaultIsZeroes;

use crate::utils::hadamard_prod;

#[derive(Clone)]
pub struct Matrix<C>(Vec<(usize, usize, C::ScalarExt)>)
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits;

impl<C> Matrix<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn mul_vector(&self, num_rows: usize, vec: &Vec<C::ScalarExt>) -> Vec<C::ScalarExt> {
        let mut result = vec![C::ScalarExt::zero(); num_rows];
        for i in 0..self.0.len() {
            let row = self.0[i].0;
            let col = self.0[i].1;
            let val = self.0[i].2;
            result[row] += val * vec[col];
        }
        result
    }
}

#[derive(Clone)]
pub struct R1CS<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub A: Matrix<C>,
    pub B: Matrix<C>,
    pub C: Matrix<C>,
    pub witness: Vec<C::ScalarExt>,
    pub public_input: Vec<C::ScalarExt>,
    pub num_cons: usize,
    pub num_vars: usize,
    pub num_input: usize,
}

impl<C> R1CS<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn produce_synthetic_r1cs(num_cons: usize, num_vars: usize, num_input: usize) -> Self {
        //        assert_eq!(num_cons, num_vars);
        let mut public_input = Vec::with_capacity(num_input);
        let mut witness = Vec::with_capacity(num_vars);

        for i in 0..num_input {
            public_input.push(C::ScalarExt::from((i + 1) as u64));
        }

        for i in 0..num_vars {
            witness.push(C::ScalarExt::from((i + 1) as u64));
        }

        let z: Vec<C::ScalarExt> = vec![public_input.clone(), witness.clone()].concat();

        let mut A: Vec<(usize, usize, C::ScalarExt)> = vec![];
        let mut B: Vec<(usize, usize, C::ScalarExt)> = vec![];
        let mut C: Vec<(usize, usize, C::ScalarExt)> = vec![];

        for i in 0..num_cons {
            let A_col = i % num_vars;
            let B_col = (i + 1) % num_vars;
            let C_col = (i + 2) % num_vars;

            // For the i'th constraint,
            // add the value 1 at the (i % num_vars)th column of A, B.
            // Compute the corresponding C_column value so that A_i * B_i = C_i
            // we apply multiplication since the Hadamard product is computed for Az ・ Bz,

            // We only _enable_ a single variable in each constraint.
            A.push((i, A_col, C::ScalarExt::one()));
            B.push((i, B_col, C::ScalarExt::one()));
            C.push((i, C_col, (z[A_col] * z[B_col]) * z[C_col].invert().unwrap()));
        }

        Self {
            A: Matrix(A),
            B: Matrix(B),
            C: Matrix(C),
            witness,
            public_input,
            num_cons,
            num_vars,
            num_input,
        }
    }

    pub fn is_sat(&self, witness: &Vec<C::ScalarExt>, public_input: &Vec<C::ScalarExt>) -> bool {
        let mut z = Vec::with_capacity(witness.len() + public_input.len() + 1);
        z.extend(public_input);
        z.extend(witness);

        let Az = self.A.mul_vector(self.num_cons, &z);
        let Bz = self.B.mul_vector(self.num_cons, &z);
        let Cz = self.C.mul_vector(self.num_cons, &z);

        hadamard_prod::<C>(&Az, &Bz) == Cz
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use halo2curves::secp256k1::Secp256k1Affine;

    #[test]
    fn test_synthetic_r1cs() {
        let num_cons = 20;
        let num_vars = 10;
        let num_input = 5;
        type C = Secp256k1Affine;

        let r1cs = R1CS::<C>::produce_synthetic_r1cs(num_cons, num_vars, num_input);

        assert_eq!(r1cs.witness.len(), num_vars);
        assert_eq!(r1cs.public_input.len(), num_input);

        assert!(r1cs.is_sat(&r1cs.witness, &r1cs.public_input));
    }
}
