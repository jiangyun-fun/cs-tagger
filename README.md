# cs-tagger

Generate [minimap2](https://github.com/lh3/minimap2)-like CS tags for BAM files.

CS tags encode the alignment between query and reference sequences in a compact string representation. This tool computes CS tags from an existing BAM alignment and a reference FASTA, then writes them as auxiliary tags into the output BAM.

## Installation

### From source

```bash
git clone https://github.com/jiangyun-fun/cs-tagger.git
cd cs-tagger
cargo build --release
```

The binary is at `target/release/cs-tagger`.

## Usage

```bash
cs-tagger -i input.bam -o output.bam reference.fa
```

### Options

```
Usage: cs-tagger [OPTIONS] --input-bam <FILE> --output-bam <FILE> <FASTA>

Arguments:
  <FASTA>              Reference FASTA file

Options:
  -i, --input-bam <FILE>    Input BAM file (use - for stdin)
  -o, --output-bam <FILE>   Output BAM file (use - for stdout)
      --add-cs <TAG>        BAM auxiliary tag name [default: cs]
      --absolute            Use absolute positions (VCF-like anchoring)
  -h, --help                Print help
  -V, --version             Print version
```

### Examples

**Basic usage** — add CS tags with default tag name `cs`:

```bash
cs-tagger -i aligned.bam -o tagged.bam reference.fa
```

**Custom tag name** — write to the `CS` tag instead:

```bash
cs-tagger -i aligned.bam -o tagged.bam --add-cs CS reference.fa
```

**Absolute positions** — anchor each operation to its 1-based reference position:

```bash
cs-tagger -i aligned.bam -o tagged.bam --absolute reference.fa
```

**Pipe mode** — read from stdin, write to stdout:

```bash
samtools view -b input.bam | cs-tagger -i - -o - reference.fa | samtools sort -o sorted.bam
```

## CS Tag Format

The CS tag encodes alignment operations as a compact string:

| Operation | Format | Description |
|-----------|--------|-------------|
| Match | `:N` | N identical bases |
| Substitution | `*xy` | ref base `x` → query base `y` |
| Insertion | `+seq` | query bases inserted (no ref consumed) |
| Deletion | `-seq` | ref bases deleted (no query consumed) |
| Intron | `~aaNag` | splice signal with intron length N |

Example CS tag: `:10*ac+gt-tgca:5:20~gt200ag`

### Absolute Position Mode

With `--absolute`, each non-match operation is prefixed by its 1-based reference position:

```
Relative:  :10*ac+gt-tgca:5
Absolute:  11*ac11+gt11-tgca16
```

## Testing

```bash
cargo test
```

## License

MIT
