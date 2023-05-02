use crate::accumulation::hadamard::{
    HadamardAccProof, HadamardAccumulator, HadamardInstance, HadamardWitness,
};
use crate::{utils::hadamard_prod, CurveAffineExt, MultiCommitGens, PRNG};
use ff::{Field, PrimeField, PrimeFieldBits};
use halo2curves::FieldExt;
use poseidon_transcript::transcript::PoseidonTranscript;
use zeroize::DefaultIsZeroes;

pub struct HadamardAccProver<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    gens: MultiCommitGens<C>,
    transcript: PoseidonTranscript<C>,
    prng: PRNG<C>,
}

impl<C> HadamardAccProver<C>
where
    C: CurveAffineExt,
    C::ScalarExt: PrimeFieldBits,
    C::ScalarExt: PrimeField<Repr = [u8; 32]>,
    C::Base: PrimeField<Repr = [u8; 32]>,
    C::ScalarExt: DefaultIsZeroes,
{
    pub fn new(gens: MultiCommitGens<C>, transcript: PoseidonTranscript<C>) -> Self {
        Self {
            gens,
            transcript,
            prng: PRNG::new(),
        }
    }

    pub fn prove_acc(
        &mut self,
        qx: &Vec<HadamardInstance<C>>,
        qw: &Vec<HadamardWitness<C>>,
    ) -> (HadamardAccumulator<C>, HadamardAccProof<C>) {
        let n = qx.len();
        let l = qw[0].a_vec.len();

        // Absorb the accumulator instances
        for acc_inst in qx {
            self.transcript.append_point(&acc_inst.0);
            self.transcript.append_point(&acc_inst.1);
            self.transcript.append_point(&acc_inst.2);
        }

        let mu = self.transcript.squeeze(1)[0];
        let mut mu_powers = vec![];
        for i in 0..n {
            mu_powers.push(mu.pow(&[i as u64, 0, 0, 0]));
        }

        let mut t_vecs = vec![Vec::with_capacity(l); 2 * n - 1];
        for i in 0..l {
            let mut a_coeffs = vec![];
            let mut b_coeffs = vec![];
            for (j, qw_i) in qw.iter().enumerate() {
                a_coeffs.push(qw_i.a_vec[i] * &mu_powers[j]);
                b_coeffs.push(qw_i.b_vec[i]);
            }
            b_coeffs.reverse();

            let mut product_coeffs = vec![C::ScalarExt::zero(); 2 * n - 1];
            for (i, a_i) in a_coeffs.iter().enumerate() {
                for (j, b_i) in b_coeffs.iter().enumerate() {
                    product_coeffs[i + j] = *a_i * *b_i;
                }
            }

            for i in 0..(2 * n - 1) {
                t_vecs[i].push(product_coeffs[i].clone());
            }
        }

        // Commit t_vecs
        let mut comm_t_vecs_low = Vec::with_capacity(n - 1);
        let mut comm_t_vecs_high = Vec::with_capacity(n - 1);
        for (i, t_vec) in t_vecs.iter().enumerate() {
            if i == n - 1 {
                continue;
            }

            if i < n {
                let c_t_i = self.gens.commit(&t_vec, &C::ScalarExt::zero());
                comm_t_vecs_low.push(c_t_i);
            } else {
                let c_t_i = self.gens.commit(&t_vec, &C::ScalarExt::zero());
                comm_t_vecs_high.push(c_t_i);
            }
        }

        let nu = self.transcript.squeeze(1)[0];

        let mut nu_powers = vec![];
        for i in 0..qx.len() {
            nu_powers.push(nu.pow(&[i as u64, 0, 0, 0]));
        }

        // Compute commitment to a(v, u)
        let mut c1 = C::identity().to_curve();
        let mut c2 = C::identity().to_curve();
        for (i, qx_i) in qx.iter().enumerate() {
            c1 += qx_i.0 * nu_powers[i] * mu_powers[i];
            c2 += qx_i.1 * nu_powers[nu_powers.len() - i - 1];
        }

        let mut c3_1 = C::identity().to_curve();
        for i in 0..comm_t_vecs_low.len() {
            c3_1 += comm_t_vecs_low[i] * nu_powers[i];
        }
        let mut c3_2 = C::identity().to_curve();
        for (i, qx_i) in qx.iter().enumerate() {
            c3_2 += qx_i.2 * mu_powers[i];
        }
        c3_2 *= nu_powers[n - 1];

        let mut c3_3 = C::identity().to_curve();
        for i in 0..comm_t_vecs_high.len() {
            c3_3 += comm_t_vecs_high[i] * nu.pow(&[(n + i) as u64, 0, 0, 0]);
        }

        let c3 = c3_1 + c3_2 + c3_3;

        // a_1 * mu^0 * nu^0 + a_2 * mu^1 * nu^1
        let mut a = vec![C::ScalarExt::zero(); l];
        let mut b = vec![C::ScalarExt::zero(); l];

        // qw_1.w_1 * mu^0 * nu^0 + qw_1.w_2 * mu^1 * nu^1
        let mut w1 = C::ScalarExt::zero();
        let mut w2 = C::ScalarExt::zero();
        let mut w3 = C::ScalarExt::zero();

        for (i, qw_i) in qw.iter().enumerate() {
            for j in 0..qw_i.a_vec.len() {
                a[j] += qw_i.a_vec[j] * &mu_powers[i] * &nu_powers[i];
            }
            w1 += qw_i.w1 * &mu_powers[i] * &nu_powers[i];

            for j in 0..qw_i.b_vec.len() {
                b[j] += qw_i.b_vec[j] * &nu_powers[n - i - 1];
            }
            w2 += qw_i.w2 * &nu_powers[n - i - 1];

            w3 += qw_i.w3 * &mu_powers[i];
        }
        w3 *= nu_powers[n - 1];

        let mut comm_t_vecs = vec![];
        comm_t_vecs.extend_from_slice(&comm_t_vecs_low);
        comm_t_vecs.extend_from_slice(&comm_t_vecs_high);

        (
            HadamardAccumulator {
                qx: HadamardInstance(c1.into(), c2.into(), c3.into()),
                qw: HadamardWitness {
                    a_vec: a,
                    b_vec: b,
                    w1,
                    w2,
                    w3,
                },
            },
            HadamardAccProof(comm_t_vecs),
        )
    }
}
