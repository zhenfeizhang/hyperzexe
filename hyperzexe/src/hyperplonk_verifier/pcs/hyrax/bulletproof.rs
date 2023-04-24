use halo2_base::utils::CurveAffineExt;

use crate::{halo2_verifier::util::transcript::TranscriptRead, EcPoint};

#[allow(non_snake_case)]
#[derive(Clone, Debug)]
struct BulletReductionProof<C: CurveAffineExt> {
    L_vec: Vec<EcPoint<C>>,
    R_vec: Vec<EcPoint<C>>,
    u_vec: Vec<EcPoint<C>>,
}

impl<C: CurveAffineExt> BulletReductionProof<C> {
    fn read_proof<T: TranscriptRead<C>>(n: usize, transcript: &mut T) -> Self {
        let mut n = n;
        let mut L_vec = Vec::new();
        let mut R_vec = Vec::new();
        let mut u_vec = Vec::new();
        while n != 1 {
            n >>= 1;
            L_vec.push(transcript.read_ec_point().unwrap());
            R_vec.push(transcript.read_ec_point().unwrap());
            u_vec.push(transcript.squeeze_challenge());
        }
        Self {
            L_vec,
            R_vec,
            u_vec,
        }
    }

    fn verification_scalars(
        &self,
    ) -> Result<
        (
            Vec<L::LoadedEcPoint>,
            Vec<L::LoadedScalar>,
            Vec<L::LoadedScalar>,
        ),
        Error,
    > {
        let loader = u_vec[0].loader();
        let lg_n = self.L_vec.len();

        // 2. Compute 1/(u_k...u_1) and 1/u_k, ..., 1/u_1
        let mut challenges_inv = self.u_vec.clone();
        L::batch_invert(challenges_inv.iter_mut());
        let allinv = loader.product(challenges_inv.iter());

        // 3. Compute u_i^2 and (1/u_i)^2
        for i in 0..lg_n {
            challenges[i] = challenges[i].square();
            challenges_inv[i] = challenges_inv[i].square();
        }
        let challenges_sq = challenges;
        let challenges_inv_sq = challenges_inv;

        // 4. Compute s values inductively.
        let mut s = Vec::with_capacity(n);
        s.push(allinv);
        for i in 1..n {
            let lg_i = (32 - 1 - (i as u32).leading_zeros()) as usize;
            let k = 1 << lg_i;
            // The challenges are stored in "creation order" as [u_k,...,u_1],
            // so u_{lg(i)+1} = is indexed by (lg_n-1) - lg_i
            let u_lg_i_sq = challenges_sq[(lg_n - 1) - lg_i];
            s.push(s[i - k] * u_lg_i_sq);
        }

        Ok((challenges_sq, challenges_inv_sq, s))
    }

    #[allow(non_snake_case)]
    fn verify(
        &self,
        a: &[L::LoadedScalar],
        Gamma: &L::LoadedEcPoint,
        G: &[L::LoadedEcPoint],
    ) -> Result<(), Error> {
        let loader = a[0].loader();
        let (u_sq, u_inv_sq, s) = self.verification_scalars()?;

        let a_hat = loader.sum_products(a.iter().zip(s.iter));
        let G_hat = loader.variable_base_msm(s.iter().zip(G.iter()));

        let points = self
            .L_vec
            .iter()
            .chain(self.R_vec.iter())
            .chain(once(Gamma));

        let scalars = u_sq
            .iter()
            .chain(u_inv_sq.iter())
            .chain(iter::once(loader.load_one()));

        let Gamma_hat = loader.variable_base_msm(scalars.iter().zip(points.iter()));

        Ok((G_hat, Gamma_hat, a_hat))
    }
}
