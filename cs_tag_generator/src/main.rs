use clap::{Arg, Command};
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
enum CsOp {
    Identical(String),           // = operation (long form)
    IdenticalLen(u32),           // : operation (short form)
    Substitution(char, char),    // * operation (ref->query)
    Insertion(String),           // + operation
    Deletion(String),            // - operation
    Intron(String, u32, String), // ~ operation (splice_start, length, splice_end)
}

impl CsOp {
    fn to_string(&self) -> String {
        match self {
            CsOp::Identical(seq) => format!("={}", seq),
            CsOp::IdenticalLen(len) => format!(":{}", len),
            CsOp::Substitution(ref_base, query_base) => {
                format!(
                    "*{}{}",
                    ref_base.to_ascii_lowercase(),
                    query_base.to_ascii_lowercase()
                )
            }
            CsOp::Insertion(seq) => format!("+{}", seq.to_lowercase()),
            CsOp::Deletion(seq) => format!("-{}", seq.to_lowercase()),
            CsOp::Intron(start, len, end) => {
                format!("~{}{}{}", start.to_lowercase(), len, end.to_lowercase())
            }
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("cs_tag_generator")
        .version("1.0")
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

    let input_bam = matches.get_one::<String>("input-bam").unwrap();
    let output_bam = matches.get_one::<String>("output-bam").unwrap();
    let reference_file = matches.get_one::<String>("reference").unwrap();
    let cs_tag_name = matches.get_one::<String>("add-cs").unwrap();
    let absolute = matches.get_flag("absolute");

    // Validate CS tag name
    if cs_tag_name.len() != 2 {
        return Err(format!(
            "CS tag name must be exactly 2 characters, got: '{}'",
            cs_tag_name
        )
        .into());
    }

    // Read reference sequences
    let reference_seqs = read_fasta(reference_file)?;

    // Process BAM file
    process_bam(input_bam, output_bam, &reference_seqs, cs_tag_name, absolute)?;

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
        if line.starts_with('>') {
            if !current_name.is_empty() {
                sequences.insert(current_name.clone(), current_seq.clone());
            }
            current_name = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_seq.clear();
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
    // let mut bam_reader = bam::Reader::from_path(input_path)?;
    let mut bam_reader = if input_path == "-" {
        bam::Reader::from_stdin()?
    } else {
        bam::Reader::from_path(input_path)?
    };
    let header = bam::Header::from_template(bam_reader.header());

    // Get header reference before the loop to avoid borrow checker issues
    let header_view = bam_reader.header().clone();

    // let mut bam_writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;
    let mut bam_writer = if output_path == "-" {
        bam::Writer::from_stdout(&header, bam::Format::Bam)?
    } else {
        bam::Writer::from_path(output_path, &header, bam::Format::Bam)?
    };

    for result in bam_reader.records() {
        let mut record = result?;

        // Skip unmapped reads
        if record.is_unmapped() {
            bam_writer.write(&record)?;
            continue;
        }

        // Get query sequence from BAM record
        let query_seq = String::from_utf8_lossy(&record.seq().as_bytes()).to_string();

        // Get reference sequence using the cloned header
        let ref_name = std::str::from_utf8(header_view.tid2name(record.tid() as u32))?;

        if let Some(ref_seq) = reference_seqs.get(ref_name) {
            // Generate CS tag
            let cs_tag = generate_cs_tag(&record, &query_seq, ref_seq)?;
            let cs_tag = if absolute {
                cs_to_absolute(&cs_tag, record.pos())
            } else {
                cs_tag
            };

            // Add CS tag to record with custom tag name
            record.push_aux(cs_tag_name.as_bytes(), bam::record::Aux::String(&cs_tag))?;
        }

        bam_writer.write(&record)?;
    }

    Ok(())
}

/// Convert a relative CS tag to absolute-position format (VCF-like anchoring).
///
/// Position advancement rules:
/// - `:N` / `=seq` → advance pos, removed from output
/// - `*xy`          → output `pos*xy`, DON'T advance
/// - `+seq`         → output `pos+seq`, DON'T advance
/// - `-seq`         → output `pos-seq`, advance by seq.len()
/// - `~aa<N>gt`     → output `pos~...`, advance by N
fn cs_to_absolute(cs_tag: &str, pos: i64) -> String {
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
                    .unwrap_or("0")
                    .parse()
                    .unwrap_or(0);
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
                // Substitution: *xy — don't advance
                i += 1;
                let ref_base = bytes[i] as char;
                i += 1;
                let query_base = bytes[i] as char;
                i += 1;
                result.push_str(&format!("{}*{}{}", abs_pos, ref_base, query_base));
            }
            b'+' => {
                // Insertion: +seq — don't advance
                i += 1;
                let start = i;
                while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                    i += 1;
                }
                let seq = std::str::from_utf8(&bytes[start..i]).unwrap_or("");
                result.push_str(&format!("{}+{}", abs_pos, seq));
            }
            b'-' => {
                // Deletion: -seq — advance by seq length
                i += 1;
                let start = i;
                while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                    i += 1;
                }
                let seq = std::str::from_utf8(&bytes[start..i]).unwrap_or("");
                result.push_str(&format!("{}-{}", abs_pos, seq));
                abs_pos += seq.len() as u64;
            }
            b'~' => {
                // Intron: ~start_seq<len>end_seq — advance by length
                i += 1;
                let ss_start = i;
                while i < bytes.len() && bytes[i].is_ascii_lowercase() {
                    i += 1;
                }
                let splice_start = std::str::from_utf8(&bytes[ss_start..i]).unwrap_or("");
                let len_start = i;
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                let intron_len: u64 = std::str::from_utf8(&bytes[len_start..i])
                    .unwrap_or("0")
                    .parse()
                    .unwrap_or(0);
                let se_start = i;
                while i < bytes.len() && bytes[i].is_ascii_lowercase() {
                    i += 1;
                }
                let splice_end = std::str::from_utf8(&bytes[se_start..i]).unwrap_or("");
                result.push_str(&format!(
                    "{}~{}{}{}",
                    abs_pos, splice_start, intron_len, splice_end
                ));
                abs_pos += intron_len;
            }
            _ => {
                i += 1;
            }
        }
    }

    result
}

fn generate_cs_tag(
    record: &bam::Record,
    query_seq: &str,
    ref_seq: &str,
) -> Result<String, Box<dyn std::error::Error>> {
    let cigar = record.cigar();
    let mut cs_ops = Vec::new();
    let mut query_pos = 0;
    let mut ref_pos = record.pos() as usize;

    let query_bytes = query_seq.as_bytes();
    let ref_bytes = ref_seq.as_bytes();

    for cigar_op in cigar.iter() {
        match cigar_op {
            bam::record::Cigar::Match(len) | bam::record::Cigar::Equal(len) => {
                let len = *len as usize;
                let mut identical_start = query_pos;
                let mut i = 0;

                while i < len {
                    if query_pos + i >= query_bytes.len() || ref_pos + i >= ref_bytes.len() {
                        break;
                    }

                    let query_base = query_bytes[query_pos + i] as char;
                    let ref_base = ref_bytes[ref_pos + i] as char;

                    if query_base.to_ascii_uppercase() != ref_base.to_ascii_uppercase() {
                        // Add any preceding identical sequence
                        if i > identical_start - query_pos {
                            let identical_len = i - (identical_start - query_pos);
                            if identical_len > 0 {
                                cs_ops.push(CsOp::IdenticalLen(identical_len as u32));
                            }
                        }

                        // Add substitution
                        cs_ops.push(CsOp::Substitution(ref_base, query_base));
                        identical_start = query_pos + i + 1;
                    }
                    i += 1;
                }

                // Add remaining identical sequence
                let remaining_len = len - (identical_start - query_pos);
                if remaining_len > 0 {
                    cs_ops.push(CsOp::IdenticalLen(remaining_len as u32));
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
                // For introns, we need splice signals (first 2 and last 2 bases)
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
            bam::record::Cigar::HardClip(_) => {
                // Hard clips don't consume query sequence
            }
            bam::record::Cigar::Pad(_) => {
                // Padding doesn't consume sequence
            }
        }
    }

    // Convert operations to string
    let cs_string: String = cs_ops.iter().map(|op| op.to_string()).collect();
    Ok(cs_string)
}
