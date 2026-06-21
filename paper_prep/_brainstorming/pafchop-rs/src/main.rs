use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

const TOOL_NAME: &str = "pafchop-rs";

#[derive(Debug, Clone)]
struct Args {
    input: String,
    length: u64,
    overlap: u64,
    keep_tags: bool,
    comparison_id: String,
    summary: Option<PathBuf>,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: "-".to_string(),
            length: 10_000,
            overlap: 0,
            keep_tags: false,
            comparison_id: "NA".to_string(),
            summary: None,
        }
    }
}

#[derive(Debug, Default, Clone)]
struct Stats {
    raw_records: u64,
    chopped_records: u64,
    raw_query_bp: u128,
    chopped_query_bp: u128,
}

fn usage() -> &'static str {
    "usage: pafchop [--input PAF|-] [--length N] [--overlap N] [--keep-tags] [--comparison-id ID] [--summary TSV]\n\
     \n\
     Streams PAF and splits each row into <=N bp query-axis fragments.\n\
     Default length is 10000 and inherited optional PAF tags are dropped to avoid\n\
     duplicating huge CIGAR/cs tags into every fragment. Input gzip is handled by\n\
     the wrapper with pigz/gzip -dc; this binary reads plain PAF from stdin/file."
}

fn parse_size(s: &str) -> Result<u64, String> {
    let s = s.trim();
    if s.is_empty() {
        return Err("empty size".to_string());
    }
    let (num, mult) = match s.as_bytes()[s.len() - 1] as char {
        'k' | 'K' => (&s[..s.len() - 1], 1_000_u64),
        'm' | 'M' => (&s[..s.len() - 1], 1_000_000_u64),
        'g' | 'G' => (&s[..s.len() - 1], 1_000_000_000_u64),
        _ => (s, 1_u64),
    };
    let value: u64 = num
        .parse()
        .map_err(|_| format!("could not parse size: {s}"))?;
    value
        .checked_mul(mult)
        .ok_or_else(|| format!("size overflows u64: {s}"))
}

fn parse_args<I>(mut it: I) -> Result<Args, String>
where
    I: Iterator<Item = String>,
{
    let mut args = Args::default();
    let _program = it.next();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                println!("{}", usage());
                std::process::exit(0);
            }
            "-i" | "--input" => args.input = it.next().ok_or("--input requires a value")?,
            "-l" | "--length" => {
                args.length = parse_size(&it.next().ok_or("--length requires a value")?)?
            }
            "-o" | "--overlap" => {
                args.overlap = parse_size(&it.next().ok_or("--overlap requires a value")?)?
            }
            "--keep-tags" => args.keep_tags = true,
            "--comparison-id" => {
                args.comparison_id = it.next().ok_or("--comparison-id requires a value")?
            }
            "--summary" => {
                args.summary = Some(PathBuf::from(
                    it.next().ok_or("--summary requires a value")?,
                ))
            }
            other if other.starts_with('-') => return Err(format!("unknown option: {other}")),
            other => args.input = other.to_string(),
        }
    }
    if args.length == 0 {
        return Err("--length must be > 0".to_string());
    }
    if args.overlap >= args.length {
        return Err("--overlap must be smaller than --length".to_string());
    }
    Ok(args)
}

fn parse_u64(fields: &[&str], idx: usize, line_no: u64) -> Result<u64, String> {
    fields
        .get(idx)
        .ok_or_else(|| format!("line {line_no}: missing PAF field {idx}"))?
        .parse::<u64>()
        .map_err(|_| {
            format!(
                "line {line_no}: invalid integer in PAF field {idx}: {}",
                fields[idx]
            )
        })
}

fn emit_fragment<W: Write>(
    out: &mut W,
    fields: &[&str],
    keep_tags: bool,
    idx: u64,
    length: u64,
    overlap: u64,
    orig_q_start: u64,
    orig_q_end: u64,
    orig_t_start: u64,
    orig_t_end: u64,
    q_start: u64,
    q_end: u64,
    t_start: u64,
    t_end: u64,
    matches: u64,
    aln_len: u64,
) -> io::Result<()> {
    for (i, field) in fields.iter().take(12).enumerate() {
        if i > 0 {
            write!(out, "\t")?;
        }
        match i {
            2 => write!(out, "{q_start}")?,
            3 => write!(out, "{q_end}")?,
            7 => write!(out, "{t_start}")?,
            8 => write!(out, "{t_end}")?,
            9 => write!(out, "{matches}")?,
            10 => write!(out, "{aln_len}")?,
            _ => write!(out, "{field}")?,
        }
    }
    if keep_tags {
        for tag in fields.iter().skip(12) {
            write!(out, "\t{tag}")?;
        }
    }
    write!(
        out,
        "\tzp:Z:{TOOL_NAME}\tzc:i:{idx}\tzl:i:{length}\tzo:i:{overlap}\tzs:i:{orig_q_start}\tze:i:{orig_q_end}\tzts:i:{orig_t_start}\tzte:i:{orig_t_end}\n"
    )
}

fn chop_line<W: Write>(
    out: &mut W,
    fields: &[&str],
    args: &Args,
    line_no: u64,
) -> Result<(u64, u128), String> {
    if fields.len() < 12 {
        return Err(format!("line {line_no}: expected at least 12 PAF fields"));
    }

    let q_start = parse_u64(fields, 2, line_no)?;
    let q_end = parse_u64(fields, 3, line_no)?;
    let t_start = parse_u64(fields, 7, line_no)?;
    let t_end = parse_u64(fields, 8, line_no)?;
    let matches = parse_u64(fields, 9, line_no)?;
    let aln_len = parse_u64(fields, 10, line_no)?;

    let span = q_end
        .checked_sub(q_start)
        .ok_or_else(|| format!("line {line_no}: q_end < q_start"))?;
    if span == 0 {
        return Ok((0, 0));
    }

    let step = args.length - args.overlap;
    let t_delta = t_end as i128 - t_start as i128;
    let mut emitted = 0_u64;
    let mut chopped_bp = 0_u128;
    let mut offset = 0_u64;

    while offset < span {
        let frag_q_start = q_start + offset;
        let frag_q_end = q_end.min(frag_q_start + args.length);
        let frag_span = frag_q_end - frag_q_start;
        if frag_span == 0 {
            break;
        }

        let frag_t_start =
            t_start as i128 + (t_delta * offset as i128 + (span / 2) as i128) / span as i128;
        let frag_t_end = t_start as i128
            + (t_delta * (offset + frag_span) as i128 + (span / 2) as i128) / span as i128;
        let frag_aln_len =
            1_u64.max((aln_len as u128 * frag_span as u128 + (span / 2) as u128) as u64 / span);
        let frag_matches = frag_aln_len
            .min((matches as u128 * frag_span as u128 + (span / 2) as u128) as u64 / span);

        emit_fragment(
            out,
            fields,
            args.keep_tags,
            emitted,
            args.length,
            args.overlap,
            q_start,
            q_end,
            t_start,
            t_end,
            frag_q_start,
            frag_q_end,
            frag_t_start.max(0) as u64,
            frag_t_end.max(0) as u64,
            frag_matches,
            frag_aln_len,
        )
        .map_err(|e| e.to_string())?;
        emitted += 1;
        chopped_bp += frag_span as u128;
        offset += step;
    }

    Ok((emitted, chopped_bp))
}

fn process<R: BufRead, W: Write>(reader: R, out: &mut W, args: &Args) -> Result<Stats, String> {
    let mut stats = Stats::default();
    for (idx, line) in reader.lines().enumerate() {
        let line_no = idx as u64 + 1;
        let line = line.map_err(|e| format!("line {line_no}: {e}"))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let q_start = parse_u64(&fields, 2, line_no)?;
        let q_end = parse_u64(&fields, 3, line_no)?;
        let raw_bp = q_end
            .checked_sub(q_start)
            .ok_or_else(|| format!("line {line_no}: q_end < q_start"))?;
        stats.raw_records += 1;
        stats.raw_query_bp += raw_bp as u128;
        let (n, bp) = chop_line(out, &fields, args, line_no)?;
        stats.chopped_records += n;
        stats.chopped_query_bp += bp;
    }
    Ok(stats)
}

fn write_summary(args: &Args, stats: &Stats) -> Result<(), String> {
    let Some(path) = &args.summary else {
        return Ok(());
    };
    let mut writer =
        BufWriter::new(File::create(path).map_err(|e| format!("{}: {e}", path.display()))?);
    writeln!(
        writer,
        "comparison_id\ttool\tversion\tchop_length_bp\toverlap_bp\tkeep_tags\traw_records\tchopped_records\traw_query_bp\tchopped_query_bp"
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        args.comparison_id,
        TOOL_NAME,
        env!("CARGO_PKG_VERSION"),
        args.length,
        args.overlap,
        args.keep_tags,
        stats.raw_records,
        stats.chopped_records,
        stats.raw_query_bp,
        stats.chopped_query_bp
    )
    .map_err(|e| e.to_string())
}

fn run() -> Result<(), String> {
    let args = parse_args(env::args())?;
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    let stats = if args.input == "-" {
        let stdin = io::stdin();
        process(BufReader::new(stdin.lock()), &mut out, &args)?
    } else {
        let file = File::open(&args.input).map_err(|e| format!("{}: {e}", args.input))?;
        process(BufReader::new(file), &mut out, &args)?
    };
    out.flush().map_err(|e| e.to_string())?;
    write_summary(&args, &stats)
}

fn main() {
    if let Err(e) = run() {
        eprintln!("pafchop: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_case(input: &str, length: u64, keep_tags: bool) -> (String, Stats) {
        let args = Args {
            length,
            keep_tags,
            comparison_id: "test".to_string(),
            ..Args::default()
        };
        let mut out = Vec::new();
        let stats = process(BufReader::new(input.as_bytes()), &mut out, &args).unwrap();
        (String::from_utf8(out).unwrap(), stats)
    }

    #[test]
    fn splits_query_axis_into_fixed_chunks() {
        let input = "q\t100\t0\t25\t+\tt\t200\t100\t125\t25\t25\t60\tcg:Z:25=\n";
        let (out, stats) = run_case(input, 10, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("q\t100\t0\t10\t+\tt\t200\t100\t110\t10\t10\t60\t"));
        assert!(lines[1].starts_with("q\t100\t10\t20\t+\tt\t200\t110\t120\t10\t10\t60\t"));
        assert!(lines[2].starts_with("q\t100\t20\t25\t+\tt\t200\t120\t125\t5\t5\t60\t"));
        assert!(lines[0].contains("zp:Z:pafchop-rs"));
        assert!(!lines[0].contains("cg:Z:25="));
        assert_eq!(stats.raw_records, 1);
        assert_eq!(stats.chopped_records, 3);
        assert_eq!(stats.raw_query_bp, 25);
        assert_eq!(stats.chopped_query_bp, 25);
    }

    #[test]
    fn can_preserve_optional_tags_when_requested() {
        let input = "q\t100\t0\t8\t+\tt\t200\t100\t108\t8\t8\t60\tcg:Z:8=\ttp:A:P\n";
        let (out, stats) = run_case(input, 10, true);
        assert_eq!(stats.chopped_records, 1);
        assert!(out.contains("cg:Z:8="));
        assert!(out.contains("tp:A:P"));
    }

    #[test]
    fn supports_overlap() {
        let args = Args {
            length: 10,
            overlap: 2,
            comparison_id: "test".to_string(),
            ..Args::default()
        };
        let input = "q\t100\t0\t25\t+\tt\t200\t100\t125\t25\t25\t60\n";
        let mut out = Vec::new();
        let stats = process(BufReader::new(input.as_bytes()), &mut out, &args).unwrap();
        let text = String::from_utf8(out).unwrap();
        let starts: Vec<&str> = text
            .lines()
            .map(|l| l.split('\t').nth(2).unwrap())
            .collect();
        assert_eq!(starts, vec!["0", "8", "16", "24"]);
        assert_eq!(stats.chopped_records, 4);
    }

    #[test]
    fn parses_metric_sizes() {
        assert_eq!(parse_size("10k").unwrap(), 10_000);
        assert_eq!(parse_size("2M").unwrap(), 2_000_000);
    }
}
