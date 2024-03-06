use std::{
    fmt::{Display, Formatter},
    ops::Deref,
    sync::{
        atomic::{self, AtomicBool, AtomicU32},
        Arc,
    },
};

/// Small helper that allows you stop lightmap generation in any time.
#[derive(Clone, Default)]
pub struct CancellationToken(pub Arc<AtomicBool>);

impl CancellationToken {
    /// Creates new cancellation token.
    pub fn new() -> Self {
        Self::default()
    }

    /// Checks if generation was cancelled.
    pub fn is_cancelled(&self) -> bool {
        self.0.load(atomic::Ordering::SeqCst)
    }

    /// Raises cancellation flag, actual cancellation is not immediate!
    pub fn cancel(&self) {
        self.0.store(true, atomic::Ordering::SeqCst)
    }
}

/// Lightmap generation stage.
#[derive(Copy, Clone, PartialOrd, PartialEq, Ord, Eq)]
#[repr(u32)]
pub enum ProgressStage {
    /// Gathering info about lights, doing precalculations.
    LightsCaching = 0,
    /// Actual lightmap generation.
    CalculatingLight = 3,
}

impl Display for ProgressStage {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            ProgressStage::LightsCaching => {
                write!(f, "Caching Lights")
            }
            ProgressStage::CalculatingLight => {
                write!(f, "Calculating Light")
            }
        }
    }
}

/// Progress internals.
#[derive(Default)]
pub struct ProgressData {
    stage: AtomicU32,
    // Range is [0; max_iterations]
    progress: AtomicU32,
    max_iterations: AtomicU32,
}

impl ProgressData {
    /// Returns progress percentage in [0; 100] range.
    pub fn progress_percent(&self) -> u32 {
        let iterations = self.max_iterations.load(atomic::Ordering::SeqCst);
        if iterations > 0 {
            self.progress.load(atomic::Ordering::SeqCst) * 100 / iterations
        } else {
            0
        }
    }

    /// Returns current stage.
    pub fn stage(&self) -> ProgressStage {
        match self.stage.load(atomic::Ordering::SeqCst) {
            0 => ProgressStage::LightsCaching,
            3 => ProgressStage::CalculatingLight,
            _ => unreachable!(),
        }
    }

    /// Sets new stage with max iterations per stage.
    pub fn set_stage(&self, stage: ProgressStage, max_iterations: u32) {
        self.max_iterations
            .store(max_iterations, atomic::Ordering::SeqCst);
        self.progress.store(0, atomic::Ordering::SeqCst);
        self.stage.store(stage as u32, atomic::Ordering::SeqCst);
    }

    /// Advances progress.
    pub fn advance_progress(&self) {
        self.progress.fetch_add(1, atomic::Ordering::SeqCst);
    }
}

/// Small helper that allows you to track progress of lightmap generation.
#[derive(Clone, Default)]
pub struct ProgressIndicator(pub Arc<ProgressData>);

impl ProgressIndicator {
    /// Creates new progress indicator.
    pub fn new() -> Self {
        Self::default()
    }
}

impl Deref for ProgressIndicator {
    type Target = ProgressData;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
