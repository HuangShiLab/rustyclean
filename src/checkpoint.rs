use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;
use dashmap::DashMap;
use tokio::fs;
use tokio::sync::Semaphore;
use tracing::warn;

use crate::error::RustycleanError;
use crate::pipeline::Checkpoint;

pub struct CheckpointManager {
    base_dir: PathBuf,
    cache: DashMap<String, Checkpoint>,
    write_semaphore: Semaphore,
}

impl CheckpointManager {
    pub async fn new(base_dir: PathBuf) -> Result<Self> {
        fs::create_dir_all(&base_dir).await?;
        Ok(Self {
            base_dir,
            cache: DashMap::new(),
            write_semaphore: Semaphore::new(100),
        })
    }

    fn checkpoint_path(&self, sample_id: &str) -> PathBuf {
        self.base_dir.join(format!("{}.json", sample_id))
    }

    pub async fn load(&self, sample_id: &str) -> Result<Option<Checkpoint>> {
        if let Some(entry) = self.cache.get(sample_id) {
            return Ok(Some(entry.clone()));
        }

        let path = self.checkpoint_path(sample_id);
        if !path.exists() {
            return Ok(None);
        }

        let content = fs::read_to_string(&path).await?;
        let checkpoint: Checkpoint = serde_json::from_str(&content).map_err(|e| {
            RustycleanError::CheckpointCorrupted(sample_id.to_string(), e.to_string())
        })?;

        if checkpoint.version != Checkpoint::CURRENT_VERSION {
            warn!(
                sample = sample_id,
                "Checkpoint version mismatch, treating as new"
            );
            return Ok(None);
        }

        self.cache
            .insert(sample_id.to_string(), checkpoint.clone());
        Ok(Some(checkpoint))
    }

    pub async fn save(&self, checkpoint: &Checkpoint) -> Result<()> {
        let _permit = self.write_semaphore.acquire().await?;

        let path = self.checkpoint_path(&checkpoint.sample_id);
        let temp_path = path.with_extension("tmp");

        let json = serde_json::to_string_pretty(checkpoint)?;
        fs::write(&temp_path, json).await?;
        fs::rename(&temp_path, &path).await?;

        self.cache
            .insert(checkpoint.sample_id.clone(), checkpoint.clone());
        Ok(())
    }

    pub async fn load_all(&self) -> Result<HashMap<String, Checkpoint>> {
        let mut entries = fs::read_dir(&self.base_dir).await?;
        let mut results = HashMap::new();

        while let Some(entry) = entries.next_entry().await? {
            let path = entry.path();
            if path.extension().is_some_and(|e| e == "json") {
                let sample_id = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or_default()
                    .to_string();

                if let Ok(Some(checkpoint)) = self.load(&sample_id).await {
                    results.insert(sample_id, checkpoint);
                }
            }
        }

        Ok(results)
    }

    pub async fn clean_completed(&self, sample_ids: &[String]) -> Result<u64> {
        let mut cleaned = 0u64;
        for id in sample_ids {
            if let Ok(Some(cp)) = self.load(id).await {
                if cp.is_complete() {
                    let path = self.checkpoint_path(id);
                    if fs::remove_file(&path).await.is_ok() {
                        self.cache.remove(id);
                        cleaned += 1;
                    }
                }
            }
        }
        Ok(cleaned)
    }
}
