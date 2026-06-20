//! Core data contracts, mirrored exactly from ANARCI's Python structures.

/// HMM state type: match / insert / delete (Python "m" / "i" / "d").
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StateType {
    M,
    I,
    D,
}

impl StateType {
    pub fn as_char(self) -> char {
        match self {
            StateType::M => 'm',
            StateType::I => 'i',
            StateType::D => 'd',
        }
    }
}

/// One element of a `state_vector`: Python `((id, type), seq_index)`.
/// `id` is the IMGT reference state 1..=128. `si` is the 0-based index into the
/// input sequence, or `None` for a delete state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct State {
    pub id: u8,
    pub typ: StateType,
    pub si: Option<usize>,
}

pub type StateVector = Vec<State>;

/// One numbered residue: Python `((position, insertion), amino_acid)`.
/// `ins` is the insertion code, always one of the static strings from the
/// alphabet / az / za tables (`" "`, `"A"`..`"Z"`, `"AA"`..`"ZZ"`).
/// `aa` is the residue byte, or `b'-'` for a gap.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Res {
    pub num: i32,
    pub ins: &'static str,
    pub aa: u8,
}

impl Res {
    pub fn new(num: i32, ins: &'static str, aa: u8) -> Self {
        Res { num, ins, aa }
    }
    /// A gap residue `((num, ' '), '-')`.
    pub fn gap(num: i32) -> Self {
        Res {
            num,
            ins: " ",
            aa: b'-',
        }
    }
}

/// Result of numbering a single domain: `(numbering, start, end)`.
/// `residues` empty means "scheme could not be applied" (e.g. CDR too long),
/// which in ANARCI is a non-error sentinel (empty list), distinct from a raised error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Numbered {
    pub residues: Vec<Res>,
    pub start: Option<usize>,
    pub end: Option<usize>,
}

/// Errors that ANARCI raises as `AssertionError`. The message is preserved
/// verbatim so the Python layer can re-raise an identical `AssertionError`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CoreError {
    /// Maps to Python `AssertionError(msg)`.
    Assertion(String),
}

impl std::fmt::Display for CoreError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CoreError::Assertion(m) => write!(f, "{}", m),
        }
    }
}

impl std::error::Error for CoreError {}

pub type CResult<T> = Result<T, CoreError>;

/// Convenience for constructing assertion errors with `format!`-style messages.
macro_rules! assertion {
    ($($arg:tt)*) => {
        $crate::types::CoreError::Assertion(format!($($arg)*))
    };
}
pub(crate) use assertion;
