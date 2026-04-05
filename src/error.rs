use crate::pipeline::PipelineStage;

#[derive(thiserror::Error, Debug)]
pub enum RustycleanError {
    #[error("Sample {sample_id} failed at stage {stage:?}: {message}")]
    PipelineFailure {
        sample_id: String,
        stage: PipelineStage,
        message: String,
    },

    #[error("Checkpoint corrupted for sample {0}: {1}")]
    CheckpointCorrupted(String, String),

    #[error("Validation failed for sample {0}: {1}")]
    ValidationFailed(String, String),

    #[error("Tool execution failed: {0}")]
    ToolExecution(String),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Serialization error: {0}")]
    Serialization(#[from] serde_json::Error),
}
