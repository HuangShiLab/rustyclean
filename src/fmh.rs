//! FracMinHash (FMH) host-removal support.
//!
//! This module implements a human-specific FracMinHash k-mer sketch used by
//! the `fmh` host-removal backend (see `pipeline.rs`) and by the standalone
//! `build-fmh-sketch` sketch-builder binary.
//!
//! Sketch construction is two-pass and memory-efficient:
//!   1. Pass 1 streams the host genome(s) and keeps every canonical k-mer
//!      whose hash passes the FMH sampling threshold (`u64::MAX / scale`).
//!   2. Pass 2 streams the microbial reference FASTA(s); any sampled microbial
//!      k-mer whose hash is also in the host set is removed (shared k-mers are
//!      treated as non-host so that conserved microbial sequences are not
//!      stripped from metagenomes). No full microbial k-mer set is ever built.
//!
//! Hashing: canonical k-mers are hashed with a base-4 polynomial rolling hash
//! (forward and reverse-complement, tracked incrementally); the canonical
//! value `min(fwd, rev)` is then mixed to a full 64-bit value with a fibonacci
//! multiply-shift avalanche before the sampling threshold is applied. The
//! mixing step is essential: raw polynomial values are confined to `4^k`,
//! far below the `u64::MAX / scale` sampling threshold for typical k.

use std::collections::HashSet;
use std::fs::File;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};

/// Sketch file magic bytes ("FMH1").
pub const SKETCH_MAGIC: &[u8; 4] = b"FMH1";
/// Sketch file format version.
pub const SKETCH_VERSION: u32 = 1;

/// Fibonacci golden-ratio constant used for multiply-shift mixing.
const GOLDEN: u64 = 0x9E37_79B9_7F4A_7C15;

/// Fast multiply-shift (fibonacci) hasher for `u64` keys. Used for the
/// sketch's hash set so lookups on already-mixed 64-bit hashes stay cheap.
#[derive(Default)]
pub struct FmhHasher(u64);

impl Hasher for FmhHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(GOLDEN);
        }
    }

    fn write_u64(&mut self, n: u64) {
        self.0 = n.wrapping_mul(GOLDEN);
        self.0 ^= self.0 >> 29;
    }
}

pub type FmhBuildHasher = BuildHasherDefault<FmhHasher>;
pub type FmhHashSet = HashSet<u64, FmhBuildHasher>;

/// Mix a (canonical) k-mer polynomial value into a full-range 64-bit hash.
#[inline]
pub fn mix64(mut x: u64) -> u64 {
    x = x.wrapping_mul(GOLDEN);
    x ^= x >> 29;
    x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x ^= x >> 32;
    x
}

/// Returns true if a (mixed) hash passes the FracMinHash sampling threshold,
/// keeping ~1/scale of all k-mers.
// Used by the build-fmh-sketch binary; unused in the main binary's copy.
#[allow(dead_code)]
#[inline]
pub fn passes_threshold(mixed: u64, scale: u32) -> bool {
    mixed < u64::MAX / scale as u64
}

/// Incremental canonical k-mer roller.
///
/// Feed one base at a time via `push`; once `k` valid bases have been seen it
/// returns the mixed canonical hash (`mix64(min(fwd, rev))`) of the current
/// window. Ambiguous bases ('N' and anything not ACGTU) reset the window.
///
/// The polynomial values fit in a `u64` for k <= 32 (max value 4^32 - 1), so
/// the wrapping arithmetic is exact for the supported k range.
pub struct KmerRoller {
    k: u32,
    pow: u64, // 4^(k-1)
    fwd: u64,
    rev: u64,
    len: u32,
}

impl KmerRoller {
    pub fn new(k: u32) -> Self {
        assert!((1..=32).contains(&k), "k must be in 1..=32, got {}", k);
        Self {
            k,
            pow: 1u64 << (2 * (k - 1)),
            fwd: 0,
            rev: 0,
            len: 0,
        }
    }

    #[inline]
    pub fn reset(&mut self) {
        self.fwd = 0;
        self.rev = 0;
        self.len = 0;
    }

    /// Push one base. Returns the mixed canonical hash of the trailing k-mer
    /// window, or `None` while the window is not yet full / was just reset.
    #[inline]
    pub fn push(&mut self, base: u8) -> Option<u64> {
        let code = match base {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            _ => {
                self.reset();
                return None;
            }
        };
        self.fwd = self.fwd.wrapping_mul(4).wrapping_add(code);
        // Sliding the reverse-complement hash: the outgoing base is the least
        // significant digit of the old rev value, so (rev >> 2) is exact.
        self.rev = (self.rev >> 2).wrapping_add((3 - code).wrapping_mul(self.pow));
        if self.len < self.k {
            self.len += 1;
        }
        if self.len == self.k {
            Some(mix64(self.fwd.min(self.rev)))
        } else {
            None
        }
    }
}

/// A human-specific FracMinHash k-mer sketch.
pub struct FmhSketch {
    pub k: u32,
    pub scale: u32,
    pub hashes: FmhHashSet,
}

// Sketch-construction helpers below are used by the build-fmh-sketch binary;
// they are dead in the main binary's copy of this module.
#[allow(dead_code)]
impl FmhSketch {
    pub fn new(k: u32, scale: u32) -> Self {
        assert!((1..=32).contains(&k), "k must be in 1..=32, got {}", k);
        assert!(scale >= 1, "scale must be >= 1");
        Self {
            k,
            scale,
            hashes: FmhHashSet::default(),
        }
    }

    /// Insert the mixed hash of one k-mer, applying the FMH sampling filter.
    #[inline]
    pub fn insert_mixed(&mut self, mixed: u64) {
        if passes_threshold(mixed, self.scale) {
            self.hashes.insert(mixed);
        }
    }

    /// Count how many k-mers of `seq` are present in the sketch.
    pub fn classify(&self, seq: &[u8]) -> u64 {
        let mut roller = KmerRoller::new(self.k);
        let mut hits = 0u64;
        for &b in seq {
            if let Some(h) = roller.push(b) {
                if self.hashes.contains(&h) {
                    hits += 1;
                }
            }
        }
        hits
    }

    /// Serialize the sketch: magic, version, k, scale, count, then the sorted
    /// little-endian u64 hashes.
    pub fn save(&self, path: &Path) -> Result<()> {
        let mut sorted: Vec<u64> = self.hashes.iter().copied().collect();
        sorted.sort_unstable();

        let file = File::create(path)
            .with_context(|| format!("failed to create sketch file {}", path.display()))?;
        let mut writer = std::io::BufWriter::new(file);
        writer.write_all(SKETCH_MAGIC)?;
        writer.write_all(&SKETCH_VERSION.to_le_bytes())?;
        writer.write_all(&self.k.to_le_bytes())?;
        writer.write_all(&self.scale.to_le_bytes())?;
        writer.write_all(&(sorted.len() as u64).to_le_bytes())?;
        for h in &sorted {
            writer.write_all(&h.to_le_bytes())?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Load and validate a sketch file written by `save`.
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("failed to open sketch file {}", path.display()))?;
        let mut reader = BufReader::new(file);

        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != SKETCH_MAGIC {
            bail!(
                "{}: not an FMH sketch (bad magic bytes)",
                path.display()
            );
        }

        let mut buf4 = [0u8; 4];
        reader.read_exact(&mut buf4)?;
        let version = u32::from_le_bytes(buf4);
        if version != SKETCH_VERSION {
            bail!(
                "{}: unsupported sketch version {} (expected {})",
                path.display(),
                version,
                SKETCH_VERSION
            );
        }

        reader.read_exact(&mut buf4)?;
        let k = u32::from_le_bytes(buf4);
        if !(1..=32).contains(&k) {
            bail!("{}: invalid k-mer size {} in sketch", path.display(), k);
        }

        reader.read_exact(&mut buf4)?;
        let scale = u32::from_le_bytes(buf4);
        if scale == 0 {
            bail!("{}: invalid scale 0 in sketch", path.display());
        }

        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8)?;
        let count = u64::from_le_bytes(buf8);

        let mut hashes = FmhHashSet::default();
        let mut prev = None;
        let mut buf8 = [0u8; 8];
        for i in 0..count {
            if reader.read_exact(&mut buf8).is_err() {
                bail!(
                    "{}: truncated sketch (expected {} hashes, found {})",
                    path.display(),
                    count,
                    i
                );
            }
            let h = u64::from_le_bytes(buf8);
            if let Some(p) = prev {
                if h <= p {
                    bail!("{}: sketch hashes are not strictly increasing", path.display());
                }
            }
            prev = Some(h);
            hashes.insert(h);
        }

        Ok(Self { k, scale, hashes })
    }
}

/// Open a (possibly gzipped) FASTA file for streaming line iteration.
#[allow(dead_code)]
pub fn open_fasta_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    use flate2::read::MultiGzDecoder;
    let file = File::open(path)
        .with_context(|| format!("failed to open FASTA {}", path.display()))?;
    if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Stream a (possibly gzipped, multi-record) FASTA file line by line.
///
/// `handle` is invoked for every sequence line; `record_start` is true when
/// the previous line was a header ('>'), i.e. at the start of a new record.
/// Header lines themselves are not passed to the handler. Callers keep a
/// `KmerRoller` and reset it whenever `record_start` is true.
#[allow(dead_code)]
pub fn stream_fasta_lines<F: FnMut(&[u8], bool)>(path: &Path, mut handle: F) -> Result<()> {
    let mut reader = open_fasta_reader(path)?;
    let mut record_start = true; // tolerate a file without a leading header
    let mut line = Vec::new();
    loop {
        line.clear();
        let n = reader
            .read_until(b'\n', &mut line)
            .with_context(|| format!("read error in FASTA {}", path.display()))?;
        if n == 0 {
            break;
        }
        let trimmed = {
            let mut t: &[u8] = &line;
            while let Some((&last, rest)) = t.split_last() {
                if last == b'\n' || last == b'\r' || last == b' ' || last == b'\t' {
                    t = rest;
                } else {
                    break;
                }
            }
            t
        };
        if trimmed.is_empty() {
            continue;
        }
        if trimmed[0] == b'>' {
            record_start = true;
            continue;
        }
        handle(trimmed, record_start);
        record_start = false;
    }
    Ok(())
}
