use clap::{Arg, Command};
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::fmt::{self, Write};
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
enum CsOp {
    IdenticalLen(u32),
    Substitution(char, char),
    Insertion(String),
    Deletion(String),
    Intron(String, u32, String),
}

impl fmt::Display for CsOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CsOp::IdenticalLen(len) => write!(f, ":{len}"),
            CsOp::Substitution(ref_base, query_base) => write!(
                f,
                "*{}{}",
                ref_base.to_ascii_lowercase(),
                query_base.to_ascii_lowercase()
            ),
            CsOp::Insertion(seq) => write!(f, "+{}", seq.to_lowercase()),
            CsOp::Deletion(seq) => write!(f, "-{}", seq.to_lowercase()),
            CsOp::Intron(start, len, end) => {
                write!(f, "~{}{}{}", start.to_lowercase(), len, end.to_lowercase())
            }
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new(env!("CARGO_PKG_NAME"))
        .version(env!("CARGO_PKG_VERSION"))
        .about("Generate minimap2-like CS tags for BAM files")
        .arg(
            Arg::new("input-bam")
                .short('i')
                .long("input-bam")
                .value_name("FILE")
                .help("Input BAM file")
                .required(true),
        )
        .arg(
            Arg::new("output-bam")
                .short('o')
                .long("output-bam")
                .value_name("FILE")
                .help("Output BAM file")
                .required(true),
        )
        .arg(
            Arg::new("reference")
                .value_name("FASTA")
                .help("Reference FASTA file")
                .required(true),
        )
        .arg(
            Arg::new("add-cs")
                .long("add-cs")
                .value_name("TAG")
                .help("CS tag name (exactly 2 characters, default: 'cs')")
                .default_value("cs"),
        )
        .arg(
            Arg::new("absolute")
                .long("absolute")
                .help("Use absolute positions in CS tag (VCF-like anchoring)")
                .action(clap::ArgAction::SetTrue),
        )
        .get_matches();

    let input_bam = matches
        .get_one::<String>("input-bam")
        .expect("input-bam is required");
    let output_bam = matches
        .get_one::<String>("output-bam")
        .expect("output-bam is required");
    let reference_file = matches
        .get_one::<String>("reference")
        .expect("reference is required");
    let cs_tag_name = matches
        .get_one::<String>("add-cs")
        .expect("add-cs has default value");
    let absolute = matches.get_flag("absolute");

    validate_tag_name(cs_tag_name)?;

    let reference_seqs = read_fasta(reference_file)?;

    process_bam(
        input_bam,
        output_bam,
        &reference_seqs,
        cs_tag_name,
        absolute,
    )?;

    Ok(())
}

/// Validate BAM auxiliary tag name per SAM spec: `[A-Za-z][A-Za-z0-9]`.
fn validate_tag_name(name: &str) -> Result<(), Box<dyn std::error::Error>> {
    if name.len() != 2 {
        return Err(format!("CS tag name must be exactly 2 characters, got: '{name}'").into());
    }
    let bytes = name.as_bytes();
    if !bytes[0].is_ascii_alphabetic() || !bytes[1].is_ascii_alphanumeric() {
        return Err(format!("CS tag name must match [A-Za-z][A-Za-z0-9], got: '{name}'").into());
    }
    Ok(())
}

fn read_fasta(file_path: &str) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut sequences = HashMap::new();
    let mut current_name = String::new();
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(name) = line.strip_prefix('>') {
            if !current_name.is_empty() {
                sequences.insert(current_name, current_seq);
            }
            current_name = name.split_whitespace().next().unwrap_or("").to_string();
            current_seq = String::new();
        } else if line.starts_with(';') {
            // FASTA comment line — skip
        } else {
            current_seq.push_str(&line);
        }
    }

    if !current_name.is_empty() {
        sequences.insert(current_name, current_seq);
    }

    Ok(sequences)
}

fn process_bam(
    input_path: &str,
    output_path: &str,
    reference_seqs: &HashMap<String, String>,
    cs_tag_name: &str,
    absolute: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut bam_reader = if input_path == "-" {
        bam::Reader::from_stdin()?
    } else {
        bam::Reader::from_path(input_path)?
    };

    let header = bam::Header::from_template(bam_reader.header());
    let header_view = bam_reader.header().clone();

    let mut bam_writer = if output_path == "-" {
        bam::Writer::from_stdout(&header, bam::Format::Bam)?
    } else {
        bam::Writer::from_path(output_path, &header, bam::Format::Bam)?
    };

    let mut missing_refs: HashMap<String, u64> = HashMap::new();

    for result in bam_reader.records() {
        let mut record = result?;

        if record.is_unmapped() {
            bam_writer.write(&record)?;
            continue;
        }

        // Validate position is non-negative before casting to usize
        let pos = record.pos();
        if pos < 0 {
            eprintln!(
                "Warning: record '{}' has negative position {pos}, skipping CS tag generation",
                String::from_utf8_lossy(record.qname()),
            );
            bam_writer.write(&record)?;
            continue;
        }

        let ref_name = std::str::from_utf8(header_view.tid2name(record.tid() as u32))?;

        if let Some(ref_seq) = reference_seqs.get(ref_name) {
            let query_bytes = record.seq().as_bytes();
            let query_str = String::from_utf8_lossy(&query_bytes);

            let cs_tag = generate_cs_tag(&record, &query_str, ref_seq, pos as usize)?;
            let cs_tag = if absolute {
                cs_to_absolute(&cs_tag, pos)?
            } else {
                cs_tag
            };

            record.push_aux(cs_tag_name.as_bytes(), bam::record::Aux::String(&cs_tag))?;
        } else {
            let count = missing_refs.entry(ref_name.to_string()).or_insert(0);
            *count += 1;
        }

        bam_writer.write(&record)?;
    }

    if !missing_refs.is_empty() {
        eprintln!("Warning: the following references were not found in the FASTA file:");
        for (name, count) in &missing_refs {
            eprintln!("  {name}: {count} records without CS tag");
        }
    }

    Ok(())
}

/// Convert a relative CS tag to absolute-position format (VCF-like anchoring).
///
/// Position advancement rules:
/// - `:N` / `=seq` → advance pos, removed from output
/// - `*xy`          → output `pos*xy`, advance by 1
/// - `+seq`         → output `pos+seq`, DON'T advance
/// - `-seq`         → output `pos-seq`, advance by seq.len()
/// - `~aa<N>gt`     → output `pos~...`, advance by N
fn cs_to_absolute(cs_tag: &str, pos: i64) -> Result<String, Box<dyn std::error::Error>> {
    let mut result = String::new();
    let mut abs_pos = (pos + 1) as u64; // 0-based BAM → 1-based
    let bytes = cs_tag.as_bytes();
    let mut i = 0;

    while i < bytes.len() {
        match bytes[i] {
            b':' => {
                // Short identical: :N
                i += 1;
                let start = i;
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                let len: u64 = std::str::from_utf8(&bytes[start..i])
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| format!("malformed :N at byte {i} in CS tag"))?;
                abs_pos += len;
            }
            b'=' => {
                // Long identical: =ACGT
                i += 1;
                let start = i;
                while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                    i += 1;
                }
                abs_pos += (i - start) as u64;
            }
            b'*' => {
                // Substitution: *xy — advance by 1 (consumes 1 ref base)
                i += 1;
                if i + 2 > bytes.len() {
                    return Err(format!("malformed *xy at byte {} in CS tag", i - 1).into());
                }
                let ref_base = bytes[i] as char;
                i += 1;
                let query_base = bytes[i] as char;
                i += 1;
                write!(result, "{abs_pos}*{ref_base}{query_base}")?;
                abs_pos += 1;
            }
            b'+' => {
                // Insertion: +seq — don't advance
                i += 1;
                let start = i;
                while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                    i += 1;
                }
                let seq = std::str::from_utf8(&bytes[start..i])
                    .ok()
                    .ok_or_else(|| format!("malformed +seq at byte {start} in CS tag"))?;
                write!(result, "{abs_pos}+{seq}")?;
            }
            b'-' => {
                // Deletion: -seq — advance by seq length
                i += 1;
                let start = i;
                while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                    i += 1;
                }
                let seq = std::str::from_utf8(&bytes[start..i])
                    .ok()
                    .ok_or_else(|| format!("malformed -seq at byte {start} in CS tag"))?;
                write!(result, "{abs_pos}-{seq}")?;
                abs_pos += seq.len() as u64;
            }
            b'~' => {
                // Intron: ~start_seq<len>end_seq — advance by length
                i += 1;
                let ss_start = i;
                while i < bytes.len() && bytes[i].is_ascii_lowercase() {
                    i += 1;
                }
                let splice_start = std::str::from_utf8(&bytes[ss_start..i])
                    .ok()
                    .ok_or_else(|| format!("malformed ~intron at byte {ss_start} in CS tag"))?;
                let len_start = i;
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                let intron_len: u64 = std::str::from_utf8(&bytes[len_start..i])
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| {
                        format!("malformed intron length at byte {len_start} in CS tag")
                    })?;
                let se_start = i;
                while i < bytes.len() && bytes[i].is_ascii_lowercase() {
                    i += 1;
                }
                let splice_end = std::str::from_utf8(&bytes[se_start..i])
                    .ok()
                    .ok_or_else(|| format!("malformed intron end at byte {se_start} in CS tag"))?;
                write!(result, "{abs_pos}~{splice_start}{intron_len}{splice_end}")?;
                abs_pos += intron_len;
            }
            _ => {
                i += 1;
            }
        }
    }

    Ok(result)
}

fn generate_cs_tag(
    record: &bam::Record,
    query_seq: &str,
    ref_seq: &str,
    mut ref_pos: usize,
) -> Result<String, Box<dyn std::error::Error>> {
    let cigar = record.cigar();
    let mut cs_ops: Vec<CsOp> = Vec::new();
    let mut query_pos = 0;

    let query_bytes = query_seq.as_bytes();
    let ref_bytes = ref_seq.as_bytes();

    for cigar_op in cigar.iter() {
        match cigar_op {
            bam::record::Cigar::Match(len) | bam::record::Cigar::Equal(len) => {
                let len = *len as usize;
                let mut identical_len: u32 = 0;

                for i in 0..len {
                    if query_pos + i >= query_bytes.len() || ref_pos + i >= ref_bytes.len() {
                        break;
                    }

                    let query_base = query_bytes[query_pos + i];
                    let ref_base = ref_bytes[ref_pos + i];

                    if !query_base.eq_ignore_ascii_case(&ref_base) {
                        if identical_len > 0 {
                            cs_ops.push(CsOp::IdenticalLen(identical_len));
                            identical_len = 0;
                        }
                        cs_ops.push(CsOp::Substitution(ref_base as char, query_base as char));
                    } else {
                        identical_len += 1;
                    }
                }

                if identical_len > 0 {
                    cs_ops.push(CsOp::IdenticalLen(identical_len));
                }

                query_pos += len;
                ref_pos += len;
            }
            bam::record::Cigar::Diff(len) => {
                let len = *len as usize;
                for i in 0..len {
                    if query_pos + i >= query_bytes.len() || ref_pos + i >= ref_bytes.len() {
                        break;
                    }

                    let query_base = query_bytes[query_pos + i] as char;
                    let ref_base = ref_bytes[ref_pos + i] as char;
                    cs_ops.push(CsOp::Substitution(ref_base, query_base));
                }
                query_pos += len;
                ref_pos += len;
            }
            bam::record::Cigar::Ins(len) => {
                let len = *len as usize;
                let end_pos = std::cmp::min(query_pos + len, query_bytes.len());
                let inserted_seq =
                    String::from_utf8_lossy(&query_bytes[query_pos..end_pos]).to_string();
                cs_ops.push(CsOp::Insertion(inserted_seq));
                query_pos += len;
            }
            bam::record::Cigar::Del(len) => {
                let len = *len as usize;
                let end_pos = std::cmp::min(ref_pos + len, ref_bytes.len());
                let deleted_seq = String::from_utf8_lossy(&ref_bytes[ref_pos..end_pos]).to_string();
                cs_ops.push(CsOp::Deletion(deleted_seq));
                ref_pos += len;
            }
            bam::record::Cigar::RefSkip(len) => {
                let len = *len as usize;
                let start_pos = ref_pos;
                let end_pos = std::cmp::min(ref_pos + len, ref_bytes.len());

                if len >= 4 && start_pos + 2 <= ref_bytes.len() && end_pos >= 2 {
                    let splice_start =
                        String::from_utf8_lossy(&ref_bytes[start_pos..start_pos + 2]).to_string();
                    let splice_end =
                        String::from_utf8_lossy(&ref_bytes[end_pos - 2..end_pos]).to_string();
                    cs_ops.push(CsOp::Intron(splice_start, len as u32, splice_end));
                }

                ref_pos += len;
            }
            bam::record::Cigar::SoftClip(len) => {
                query_pos += *len as usize;
            }
            bam::record::Cigar::HardClip(_) => {}
            bam::record::Cigar::Pad(_) => {}
        }
    }

    let cs_string: String = cs_ops.iter().map(|op| op.to_string()).collect();
    Ok(cs_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- validate_tag_name ---

    #[test]
    fn valid_tag_names() {
        assert!(validate_tag_name("cs").is_ok());
        assert!(validate_tag_name("CS").is_ok());
        assert!(validate_tag_name("Z0").is_ok());
        assert!(validate_tag_name("nm").is_ok());
    }

    #[test]
    fn invalid_tag_name_wrong_length() {
        assert!(validate_tag_name("c").is_err());
        assert!(validate_tag_name("cst").is_err());
        assert!(validate_tag_name("").is_err());
    }

    #[test]
    fn invalid_tag_name_wrong_chars() {
        assert!(validate_tag_name("0c").is_err());
        assert!(validate_tag_name("!c").is_err());
        assert!(validate_tag_name("c!").is_err());
    }

    // --- CsOp Display ---

    #[test]
    fn cs_op_display() {
        assert_eq!(CsOp::IdenticalLen(5).to_string(), ":5");
        assert_eq!(CsOp::Substitution('A', 'G').to_string(), "*ag");
        assert_eq!(CsOp::Insertion("ACGT".into()).to_string(), "+acgt");
        assert_eq!(CsOp::Deletion("TGCA".into()).to_string(), "-tgca");
        assert_eq!(
            CsOp::Intron("gt".into(), 100, "ag".into()).to_string(),
            "~gt100ag"
        );
    }

    // --- cs_to_absolute ---

    #[test]
    fn cs_to_absolute_substitution() {
        let result = cs_to_absolute("*ac", 10).unwrap();
        assert_eq!(result, "11*ac");
    }

    #[test]
    fn cs_to_absolute_insertion() {
        let result = cs_to_absolute("+acg", 10).unwrap();
        assert_eq!(result, "11+acg");
    }

    #[test]
    fn cs_to_absolute_deletion() {
        let result = cs_to_absolute("-acg", 10).unwrap();
        assert_eq!(result, "11-acg");
    }

    #[test]
    fn cs_to_absolute_intron() {
        let result = cs_to_absolute("~gt100ag", 10).unwrap();
        assert_eq!(result, "11~gt100ag");
    }

    #[test]
    fn cs_to_absolute_short_identical_advances_pos() {
        let result = cs_to_absolute(":5*ac", 10).unwrap();
        assert_eq!(result, "16*ac"); // 11 + 5 = 16, * at 16
    }

    #[test]
    fn cs_to_absolute_long_identical_advances_pos() {
        let result = cs_to_absolute("=ACGT*ac", 10).unwrap();
        assert_eq!(result, "15*ac"); // 11 + 4 = 15, * at 15
    }

    #[test]
    fn cs_to_absolute_mixed_operations() {
        // :10 → advance to 110, no output
        // +ac → "110+ac", no advance (insertion doesn't consume ref)
        // -tgca → "110-tgca", advance by 4 to 114
        // :5 → advance to 119, no output
        // *gt → "119*gt", advance by 1 to 120
        // ~gt200ag → "120~gt200ag", advance by 200
        let result = cs_to_absolute(":10+ac-tgca:5*gt~gt200ag", 99).unwrap();
        assert_eq!(result, "110+ac110-tgca119*gt120~gt200ag");
    }

    #[test]
    fn cs_to_absolute_consecutive_substitutions() {
        // Each *xy consumes 1 ref base → positions must increment
        let result = cs_to_absolute("*ac*tc*ag", 52).unwrap();
        assert_eq!(result, "53*ac54*tc55*ag");
    }

    #[test]
    fn cs_to_absolute_empty() {
        let result = cs_to_absolute("", 10).unwrap();
        assert_eq!(result, "");
    }

    #[test]
    fn cs_to_absolute_malformed_star() {
        assert!(cs_to_absolute("*", 10).is_err());
        assert!(cs_to_absolute("*a", 10).is_err());
    }

    #[test]
    fn cs_to_absolute_malformed_colon_no_digits() {
        assert!(cs_to_absolute(":", 10).is_err());
    }

    // --- read_fasta ---

    fn temp_dir(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(name);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn read_fasta_multi_sequence() {
        let dir = temp_dir("cs_tag_test_multi");
        let path = dir.join("test.fasta");
        std::fs::write(&path, ">chr1\nACGT\nTGCA\n>chr2\nAAAA\n").unwrap();

        let seqs = read_fasta(path.to_str().unwrap()).unwrap();
        assert_eq!(seqs.get("chr1").unwrap(), "ACGTTGCA");
        assert_eq!(seqs.get("chr2").unwrap(), "AAAA");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn read_fasta_with_comments() {
        let dir = temp_dir("cs_tag_test_comments");
        let path = dir.join("test.fasta");
        std::fs::write(
            &path,
            ">chr1 some description\nACGT\n; comment line\nTGCA\n",
        )
        .unwrap();

        let seqs = read_fasta(path.to_str().unwrap()).unwrap();
        assert_eq!(seqs.get("chr1").unwrap(), "ACGTTGCA");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn read_fasta_single_sequence() {
        let dir = temp_dir("cs_tag_test_single");
        let path = dir.join("test.fasta");
        std::fs::write(&path, ">solo\nAAAA\n").unwrap();

        let seqs = read_fasta(path.to_str().unwrap()).unwrap();
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs.get("solo").unwrap(), "AAAA");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn read_fasta_empty_file() {
        let dir = temp_dir("cs_tag_test_empty");
        let path = dir.join("test.fasta");
        std::fs::write(&path, "").unwrap();

        let seqs = read_fasta(path.to_str().unwrap()).unwrap();
        assert!(seqs.is_empty());

        std::fs::remove_dir_all(&dir).ok();
    }
}
