use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use xxhash_rust::xxh3::xxh3_64;

#[derive(Debug, Clone)]
pub struct Sample {
    pub id: String,
    pub r1: PathBuf,
    pub r2: Option<PathBuf>,
    pub output_dir: PathBuf,
}

impl Sample {
    pub fn is_paired(&self) -> bool {
        self.r2.is_some()
    }

    /// Compute a fast hash based on file metadata (size + mtime) for checkpoint invalidation.
    pub async fn compute_hash(&self) -> Result<u64> {
        let m1 = tokio::fs::metadata(&self.r1).await?;
        let mut input = format!(
            "{}:{}:{}",
            self.r1.display(),
            m1.len(),
            m1.modified()?.duration_since(std::time::UNIX_EPOCH)?.as_secs(),
        );

        if let Some(r2) = &self.r2 {
            let m2 = tokio::fs::metadata(r2).await?;
            input.push_str(&format!(
                ":{}:{}:{}",
                r2.display(),
                m2.len(),
                m2.modified()?.duration_since(std::time::UNIX_EPOCH)?.as_secs(),
            ));
        }

        Ok(xxh3_64(input.as_bytes()))
    }
}

/// Parse a sample list file (TSV, no header required).
///
/// Supported formats:
/// - 2 columns: sample_id \t r1_path           (single-end)
/// - 3 columns: sample_id \t r1_path \t r2_path (paired-end)
///
/// Lines starting with '#' are treated as comments and skipped.
pub async fn parse_sample_list(path: &Path, output_base: &Path) -> Result<Vec<Sample>> {
    let content = tokio::fs::read_to_string(path)
        .await
        .with_context(|| format!("Failed to read sample list: {}", path.display()))?;

    let mut samples = Vec::new();

    for (line_num, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        let sample = match fields.len() {
            2 => Sample {
                id: fields[0].to_string(),
                r1: PathBuf::from(fields[1]),
                r2: None,
                output_dir: output_base.join(fields[0]),
            },
            3 => Sample {
                id: fields[0].to_string(),
                r1: PathBuf::from(fields[1]),
                r2: Some(PathBuf::from(fields[2])),
                output_dir: output_base.join(fields[0]),
            },
            _ => bail!(
                "{}:{}: expected 2 or 3 tab-separated columns, got {}",
                path.display(),
                line_num + 1,
                fields.len()
            ),
        };

        samples.push(sample);
    }

    Ok(samples)
}

/// Create a single sample from direct CLI input (--r1, optional --r2).
///
/// The sample ID is derived from the parent directory of R1.  This avoids the
/// common case where every input file is named `reads_R1.fastq.gz` and would
/// otherwise collide on the generic id "reads".
pub fn sample_from_paths(r1: PathBuf, r2: Option<PathBuf>, output_base: &Path) -> Result<Sample> {
    let id = r1
        .parent()
        .and_then(|p| p.file_name())
        .and_then(|s| s.to_str())
        .unwrap_or("sample")
        .to_string();

    Ok(Sample {
        id: id.clone(),
        r1,
        r2,
        output_dir: output_base.join(&id),
    })
}

/// Validate that all input files exist.
pub async fn validate_inputs(samples: &[Sample]) -> Result<()> {
    let mut errors = Vec::new();

    for sample in samples {
        if !sample.r1.exists() {
            errors.push(format!("{}: R1 file not found: {}", sample.id, sample.r1.display()));
        }
        if let Some(r2) = &sample.r2 {
            if !r2.exists() {
                errors.push(format!("{}: R2 file not found: {}", sample.id, r2.display()));
            }
        }
    }

    if !errors.is_empty() {
        bail!("Input validation failed:\n  {}", errors.join("\n  "));
    }

    Ok(())
}
