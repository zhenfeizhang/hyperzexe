use serde::{Serialize, Deserialize};

mod hyrax_verifier;
mod sumcheck_verifier;
mod verifier;

#[derive(Clone, Debug)]
pub enum Error {
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Protocol<C, L = loader::native::NativeLoader>
where
    C: util::arithmetic::CurveAffine,
    L: loader::Loader<C>,
{
    #[serde(bound(
        serialize = "L::LoadedEcPoint: Serialize",
        deserialize = "L::LoadedEcPoint: Deserialize<'de>"
    ))]
    pub num_instance: Vec<usize>,
    pub preprocessed: Vec<L::LoadedEcPoint>,
    pub num_witness: Vec<usize>,
    // permutations & lookup

    
    pub num_challenge: Vec<usize>,
    pub evaluations: Vec<util::protocol::Query>,
    pub queries: Vec<util::protocol::Query>,
    // Minor customization
    #[serde(bound(
        serialize = "L::LoadedScalar: Serialize",
        deserialize = "L::LoadedScalar: Deserialize<'de>"
    ))]
    pub transcript_initial_state: Option<L::LoadedScalar>,
    pub instance_committing_key: Option<util::protocol::InstanceCommittingKey<C>>,
}
