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
    chunk_mode: ChunkMode,
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
            chunk_mode: ChunkMode::RowStart,
            keep_tags: false,
            comparison_id: "NA".to_string(),
            summary: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChunkMode {
    RowStart,
    QueryGrid,
}

impl ChunkMode {
    fn as_str(self) -> &'static str {
        match self {
            ChunkMode::RowStart => "row-start",
            ChunkMode::QueryGrid => "query-grid",
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OpKind {
    Match,
    Equal,
    Diff,
    Ins,
    Del,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CigarOp {
    len: u64,
    kind: OpKind,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum CsOp {
    MatchLen(u64),
    Equal(String),
    Diff(String),
    Ins(String),
    Del(String),
}

#[derive(Debug, Default, Clone, PartialEq)]
struct ChunkMetrics {
    matches: u64,
    aln_len: u64,
    nm: u64,
    gap_compressed_nm: u64,
}

#[derive(Debug, Default, Clone)]
struct ChunkPlan {
    q_start: u64,
    q_end: u64,
    t_min: Option<u64>,
    t_max: Option<u64>,
    cigar: Vec<CigarOp>,
    cs: Option<Vec<CsOp>>,
    metrics: ChunkMetrics,
}

impl ChunkPlan {
    fn new(q_start: u64, q_end: u64, has_cs: bool) -> Self {
        Self {
            q_start,
            q_end,
            t_min: None,
            t_max: None,
            cigar: Vec::new(),
            cs: has_cs.then(Vec::new),
            metrics: ChunkMetrics::default(),
        }
    }

    fn add_target_interval(&mut self, lo: u64, hi: u64) {
        if lo == hi {
            return;
        }
        self.t_min = Some(self.t_min.map_or(lo, |x| x.min(lo)));
        self.t_max = Some(self.t_max.map_or(hi, |x| x.max(hi)));
    }

    fn push_cigar(&mut self, op: CigarOp) {
        push_cigar_op(&mut self.cigar, op);
    }

    fn push_cs(&mut self, op: CsOp) {
        if let Some(cs) = &mut self.cs {
            push_cs_op(cs, op);
        }
    }
}

fn usage() -> &'static str {
    "usage: pafchop [--input PAF|-] [--length N] [--overlap N] [--chunk-mode row-start|query-grid] [--query-grid] [--keep-tags] [--comparison-id ID] [--summary TSV]\n\
     \n\
     Streams PAF and splits each row into <=N bp query-axis fragments.\n\
     Default chunk mode is row-start, which starts chunks at each row's q_start.\n\
     query-grid mode uses absolute query-coordinate boundaries at 0,N,2N,...\n\
     and currently requires --overlap 0.\n\
     Exact chunking requires cg:Z CIGAR operations. Alignment-derived columns\n\
     and tags are recomputed per fragment; stale cg/cs/NM/dv/de/df tags are\n\
     never copied from the input row."
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
            "--chunk-mode" => {
                let value = it.next().ok_or("--chunk-mode requires a value")?;
                args.chunk_mode = match value.as_str() {
                    "row-start" | "row_start" | "row" => ChunkMode::RowStart,
                    "query-grid" | "query_grid" | "query" => ChunkMode::QueryGrid,
                    _ => {
                        return Err(format!(
                            "--chunk-mode must be row-start or query-grid, got {value}"
                        ))
                    }
                };
            }
            "--query-grid" => args.chunk_mode = ChunkMode::QueryGrid,
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
    if args.chunk_mode == ChunkMode::QueryGrid && args.overlap != 0 {
        return Err("--chunk-mode query-grid currently requires --overlap 0".to_string());
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

fn parse_strand(s: &str, line_no: u64) -> Result<Strand, String> {
    match s {
        "+" => Ok(Strand::Forward),
        "-" => Ok(Strand::Reverse),
        _ => Err(format!("line {line_no}: invalid strand: {s}")),
    }
}

fn q_consumes(kind: OpKind) -> bool {
    matches!(
        kind,
        OpKind::Match | OpKind::Equal | OpKind::Diff | OpKind::Ins
    )
}

fn t_consumes(kind: OpKind) -> bool {
    matches!(
        kind,
        OpKind::Match | OpKind::Equal | OpKind::Diff | OpKind::Del
    )
}

fn op_code(kind: OpKind) -> char {
    match kind {
        OpKind::Match => 'M',
        OpKind::Equal => '=',
        OpKind::Diff => 'X',
        OpKind::Ins => 'I',
        OpKind::Del => 'D',
    }
}

fn parse_cigar(cigar: &str, line_no: u64) -> Result<Vec<CigarOp>, String> {
    let mut ops = Vec::new();
    let mut len = 0_u64;
    let mut saw_digit = false;
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            saw_digit = true;
            len = len
                .checked_mul(10)
                .and_then(|x| x.checked_add(ch as u64 - '0' as u64))
                .ok_or_else(|| format!("line {line_no}: CIGAR length overflow"))?;
            continue;
        }
        if !saw_digit || len == 0 {
            return Err(format!("line {line_no}: malformed cg:Z CIGAR"));
        }
        let kind = match ch {
            'M' => OpKind::Match,
            '=' => OpKind::Equal,
            'X' => OpKind::Diff,
            'I' => OpKind::Ins,
            'D' => OpKind::Del,
            _ => return Err(format!("line {line_no}: unsupported CIGAR op {ch}")),
        };
        ops.push(CigarOp { len, kind });
        len = 0;
        saw_digit = false;
    }
    if saw_digit {
        return Err(format!("line {line_no}: trailing CIGAR length without op"));
    }
    if ops.is_empty() {
        return Err(format!("line {line_no}: empty cg:Z CIGAR"));
    }
    Ok(ops)
}

fn parse_cs(cs: &str, line_no: u64) -> Result<Vec<CsOp>, String> {
    let bytes = cs.as_bytes();
    let mut idx = 0_usize;
    let mut ops = Vec::new();
    while idx < bytes.len() {
        match bytes[idx] as char {
            ':' => {
                idx += 1;
                let start = idx;
                while idx < bytes.len() && bytes[idx].is_ascii_digit() {
                    idx += 1;
                }
                if start == idx {
                    return Err(format!("line {line_no}: malformed cs match length"));
                }
                let len = cs[start..idx]
                    .parse::<u64>()
                    .map_err(|_| format!("line {line_no}: cs match length overflow"))?;
                ops.push(CsOp::MatchLen(len));
            }
            '=' | '+' | '-' => {
                let op = bytes[idx] as char;
                idx += 1;
                let start = idx;
                while idx < bytes.len()
                    && !matches!(bytes[idx] as char, ':' | '=' | '*' | '+' | '-' | '~')
                {
                    idx += 1;
                }
                if start == idx {
                    return Err(format!("line {line_no}: empty cs operation"));
                }
                let seq = cs[start..idx].to_string();
                match op {
                    '=' => ops.push(CsOp::Equal(seq)),
                    '+' => ops.push(CsOp::Ins(seq)),
                    '-' => ops.push(CsOp::Del(seq)),
                    _ => unreachable!(),
                }
            }
            '*' => {
                idx += 1;
                if idx + 2 > bytes.len() {
                    return Err(format!("line {line_no}: short cs mismatch"));
                }
                let seq = cs[idx..idx + 2].to_string();
                idx += 2;
                ops.push(CsOp::Diff(seq));
            }
            '~' => return Err(format!("line {line_no}: unsupported cs intron operation")),
            other => return Err(format!("line {line_no}: malformed cs operation {other}")),
        }
    }
    Ok(ops)
}

fn find_tag<'a>(fields: &'a [&str], prefix: &str) -> Option<&'a str> {
    fields
        .iter()
        .skip(12)
        .find_map(|tag| tag.strip_prefix(prefix))
}

fn has_tag(fields: &[&str], prefix: &str) -> bool {
    fields.iter().skip(12).any(|tag| tag.starts_with(prefix))
}

fn push_cigar_op(ops: &mut Vec<CigarOp>, op: CigarOp) {
    if op.len == 0 {
        return;
    }
    if let Some(last) = ops.last_mut() {
        if last.kind == op.kind {
            last.len += op.len;
            return;
        }
    }
    ops.push(op);
}

fn push_cs_op(ops: &mut Vec<CsOp>, op: CsOp) {
    match op {
        CsOp::MatchLen(0) => {}
        CsOp::Equal(seq) | CsOp::Ins(seq) | CsOp::Del(seq) if seq.is_empty() => {}
        other => ops.push(other),
    }
}

fn cigar_to_string(ops: &[CigarOp]) -> String {
    let mut out = String::new();
    for op in ops {
        out.push_str(&op.len.to_string());
        out.push(op_code(op.kind));
    }
    out
}

fn cs_to_string(ops: &[CsOp]) -> String {
    let mut out = String::new();
    for op in ops {
        match op {
            CsOp::MatchLen(len) => {
                out.push(':');
                out.push_str(&len.to_string());
            }
            CsOp::Equal(seq) => {
                out.push('=');
                out.push_str(seq);
            }
            CsOp::Diff(seq) => {
                out.push('*');
                out.push_str(seq);
            }
            CsOp::Ins(seq) => {
                out.push('+');
                out.push_str(seq);
            }
            CsOp::Del(seq) => {
                out.push('-');
                out.push_str(seq);
            }
        }
    }
    out
}

fn target_interval(pos: u64, len: u64, strand: Strand) -> (u64, u64, u64) {
    match strand {
        Strand::Forward => (pos, pos + len, pos + len),
        Strand::Reverse => (pos - len, pos, pos - len),
    }
}

fn add_metrics(metrics: &mut ChunkMetrics, kind: OpKind, len: u64) {
    metrics.aln_len += len;
    match kind {
        OpKind::Match | OpKind::Equal => metrics.matches += len,
        OpKind::Diff => {
            metrics.nm += len;
            metrics.gap_compressed_nm += len;
        }
        OpKind::Ins | OpKind::Del => {
            metrics.nm += len;
            metrics.gap_compressed_nm += 1;
        }
    }
}

fn metrics_from_cs(ops: &[CsOp]) -> ChunkMetrics {
    let mut metrics = ChunkMetrics::default();
    for op in ops {
        match op {
            CsOp::MatchLen(len) => {
                metrics.matches += *len;
                metrics.aln_len += *len;
            }
            CsOp::Equal(seq) => {
                metrics.matches += seq.len() as u64;
                metrics.aln_len += seq.len() as u64;
            }
            CsOp::Diff(_) => {
                metrics.aln_len += 1;
                metrics.nm += 1;
                metrics.gap_compressed_nm += 1;
            }
            CsOp::Ins(seq) => {
                metrics.aln_len += seq.len() as u64;
                metrics.nm += seq.len() as u64;
                metrics.gap_compressed_nm += 1;
            }
            CsOp::Del(seq) => {
                metrics.aln_len += seq.len() as u64;
                metrics.nm += seq.len() as u64;
                metrics.gap_compressed_nm += 1;
            }
        }
    }
    metrics
}

fn op_chunk_for_deletion(
    q_pos: u64,
    q_start: u64,
    q_end: u64,
    chunks: &[ChunkPlan],
) -> Option<usize> {
    chunks.iter().position(|chunk| {
        (q_pos >= chunk.q_start && q_pos < chunk.q_end)
            || (q_pos == q_end && q_pos == chunk.q_end)
            || (q_pos == q_start && q_pos == chunk.q_start)
    })
}

fn slice_ascii(s: &str, start: u64, len: u64) -> String {
    s.as_bytes()[start as usize..(start + len) as usize]
        .iter()
        .map(|b| *b as char)
        .collect()
}

fn apply_cs_to_chunks(
    cs_ops: &[CsOp],
    q_start: u64,
    q_end: u64,
    chunks: &mut [ChunkPlan],
) -> Result<(), String> {
    let mut q_pos = q_start;
    for op in cs_ops {
        match op {
            CsOp::MatchLen(len) => {
                let op_q_start = q_pos;
                let op_q_end = q_pos + *len;
                for chunk in chunks.iter_mut() {
                    let ov_start = op_q_start.max(chunk.q_start);
                    let ov_end = op_q_end.min(chunk.q_end);
                    if ov_start < ov_end {
                        chunk.push_cs(CsOp::MatchLen(ov_end - ov_start));
                    }
                }
                q_pos = op_q_end;
            }
            CsOp::Equal(seq) => {
                let len = seq.len() as u64;
                let op_q_start = q_pos;
                let op_q_end = q_pos + len;
                for chunk in chunks.iter_mut() {
                    let ov_start = op_q_start.max(chunk.q_start);
                    let ov_end = op_q_end.min(chunk.q_end);
                    if ov_start < ov_end {
                        chunk.push_cs(CsOp::Equal(slice_ascii(
                            seq,
                            ov_start - op_q_start,
                            ov_end - ov_start,
                        )));
                    }
                }
                q_pos = op_q_end;
            }
            CsOp::Diff(seq) => {
                let op_q_start = q_pos;
                let op_q_end = q_pos + 1;
                for chunk in chunks.iter_mut() {
                    if op_q_start >= chunk.q_start && op_q_start < chunk.q_end {
                        chunk.push_cs(CsOp::Diff(seq.clone()));
                    }
                }
                q_pos = op_q_end;
            }
            CsOp::Ins(seq) => {
                let len = seq.len() as u64;
                let op_q_start = q_pos;
                let op_q_end = q_pos + len;
                for chunk in chunks.iter_mut() {
                    let ov_start = op_q_start.max(chunk.q_start);
                    let ov_end = op_q_end.min(chunk.q_end);
                    if ov_start < ov_end {
                        chunk.push_cs(CsOp::Ins(slice_ascii(
                            seq,
                            ov_start - op_q_start,
                            ov_end - ov_start,
                        )));
                    }
                }
                q_pos = op_q_end;
            }
            CsOp::Del(seq) => {
                if let Some(idx) = op_chunk_for_deletion(q_pos, q_start, q_end, chunks) {
                    chunks[idx].push_cs(CsOp::Del(seq.clone()));
                }
            }
        }
    }
    if q_pos != q_end {
        return Err(format!(
            "cs:Z query span {q_pos} does not match PAF query end {q_end}"
        ));
    }
    Ok(())
}

fn plan_chunks(fields: &[&str], args: &Args, line_no: u64) -> Result<Vec<ChunkPlan>, String> {
    let q_start = parse_u64(fields, 2, line_no)?;
    let q_end = parse_u64(fields, 3, line_no)?;
    let strand = parse_strand(
        fields
            .get(4)
            .ok_or_else(|| format!("line {line_no}: missing strand"))?,
        line_no,
    )?;
    let t_start = parse_u64(fields, 7, line_no)?;
    let t_end = parse_u64(fields, 8, line_no)?;
    let cg = find_tag(fields, "cg:Z:")
        .ok_or_else(|| format!("line {line_no}: cannot chop exactly without cg:Z CIGAR"))?;
    let cigar = parse_cigar(cg, line_no)?;
    let cs_ops = find_tag(fields, "cs:Z:")
        .map(|cs| parse_cs(cs, line_no))
        .transpose()?;
    if cs_ops.is_none() && cigar.iter().any(|op| op.kind == OpKind::Match) {
        return Err(format!(
            "line {line_no}: cg:Z uses ambiguous M operations; cs:Z is required to recompute identity exactly"
        ));
    }

    let span = q_end
        .checked_sub(q_start)
        .ok_or_else(|| format!("line {line_no}: q_end < q_start"))?;
    if span == 0 {
        return Ok(Vec::new());
    }

    let mut chunks = Vec::new();
    match args.chunk_mode {
        ChunkMode::RowStart => {
            let step = args.length - args.overlap;
            let mut offset = 0_u64;
            while offset < span {
                let chunk_q_start = q_start + offset;
                let chunk_q_end = q_end.min(chunk_q_start + args.length);
                chunks.push(ChunkPlan::new(chunk_q_start, chunk_q_end, cs_ops.is_some()));
                offset += step;
            }
        }
        ChunkMode::QueryGrid => {
            let mut grid_start = (q_start / args.length) * args.length;
            while grid_start < q_end {
                let grid_end = grid_start.saturating_add(args.length);
                let chunk_q_start = q_start.max(grid_start);
                let chunk_q_end = q_end.min(grid_end);
                if chunk_q_start < chunk_q_end {
                    chunks.push(ChunkPlan::new(chunk_q_start, chunk_q_end, cs_ops.is_some()));
                }
                if grid_end == u64::MAX {
                    break;
                }
                grid_start = grid_end;
            }
        }
    }

    let mut q_pos = q_start;
    let mut t_pos = match strand {
        Strand::Forward => t_start,
        Strand::Reverse => t_end,
    };
    for op in &cigar {
        if q_consumes(op.kind) {
            let op_q_start = q_pos;
            let op_q_end = q_pos + op.len;
            for chunk in chunks.iter_mut() {
                let ov_start = op_q_start.max(chunk.q_start);
                let ov_end = op_q_end.min(chunk.q_end);
                if ov_start < ov_end {
                    let ov_len = ov_end - ov_start;
                    chunk.push_cigar(CigarOp {
                        len: ov_len,
                        kind: op.kind,
                    });
                    add_metrics(&mut chunk.metrics, op.kind, ov_len);
                    if t_consumes(op.kind) {
                        let t_offset = ov_start - op_q_start;
                        let oriented_start = match strand {
                            Strand::Forward => t_pos + t_offset,
                            Strand::Reverse => t_pos - t_offset,
                        };
                        let (lo, hi, _) = target_interval(oriented_start, ov_len, strand);
                        chunk.add_target_interval(lo, hi);
                    }
                }
            }
            q_pos = op_q_end;
        } else if op.kind == OpKind::Del {
            let idx = op_chunk_for_deletion(q_pos, q_start, q_end, &chunks)
                .ok_or_else(|| format!("line {line_no}: deletion outside query span"))?;
            chunks[idx].push_cigar(*op);
            add_metrics(&mut chunks[idx].metrics, op.kind, op.len);
            let (lo, hi, _) = target_interval(t_pos, op.len, strand);
            chunks[idx].add_target_interval(lo, hi);
        }

        if t_consumes(op.kind) {
            let (_, _, next) = target_interval(t_pos, op.len, strand);
            t_pos = next;
        }
    }
    if q_pos != q_end {
        return Err(format!(
            "line {line_no}: cg:Z query span ended at {q_pos}, expected {q_end}"
        ));
    }
    match strand {
        Strand::Forward if t_pos != t_end => {
            return Err(format!(
                "line {line_no}: cg:Z target span ended at {t_pos}, expected {t_end}"
            ));
        }
        Strand::Reverse if t_pos != t_start => {
            return Err(format!(
                "line {line_no}: cg:Z target span ended at {t_pos}, expected {t_start}"
            ));
        }
        _ => {}
    }

    if let Some(cs_ops) = &cs_ops {
        apply_cs_to_chunks(cs_ops, q_start, q_end, &mut chunks)
            .map_err(|e| format!("line {line_no}: {e}"))?;
        for chunk in &mut chunks {
            chunk.metrics = metrics_from_cs(chunk.cs.as_deref().unwrap_or(&[]));
        }
    }

    Ok(chunks)
}

fn emit_fragment<W: Write>(
    out: &mut W,
    fields: &[&str],
    args: &Args,
    idx: u64,
    orig_q_start: u64,
    orig_q_end: u64,
    orig_t_start: u64,
    orig_t_end: u64,
    chunk: &ChunkPlan,
) -> io::Result<()> {
    let t_start = chunk.t_min.unwrap_or(orig_t_start);
    let t_end = chunk.t_max.unwrap_or(orig_t_start);
    for (i, field) in fields.iter().take(12).enumerate() {
        if i > 0 {
            write!(out, "\t")?;
        }
        match i {
            2 => write!(out, "{}", chunk.q_start)?,
            3 => write!(out, "{}", chunk.q_end)?,
            7 => write!(out, "{t_start}")?,
            8 => write!(out, "{t_end}")?,
            9 => write!(out, "{}", chunk.metrics.matches)?,
            10 => write!(out, "{}", chunk.metrics.aln_len)?,
            _ => write!(out, "{field}")?,
        }
    }

    if args.keep_tags {
        for tag in fields.iter().skip(12) {
            if is_alignment_derived_tag(tag) {
                continue;
            }
            write!(out, "\t{tag}")?;
        }
    }

    write!(out, "\tcg:Z:{}", cigar_to_string(&chunk.cigar))?;
    if let Some(cs) = &chunk.cs {
        write!(out, "\tcs:Z:{}", cs_to_string(cs))?;
    }
    if has_tag(fields, "NM:i:") {
        write!(out, "\tNM:i:{}", chunk.metrics.nm)?;
    }
    if has_tag(fields, "dv:f:") {
        write!(
            out,
            "\tdv:f:{:.6}",
            divergence(chunk.metrics.nm, chunk.metrics.aln_len)
        )?;
    }
    if has_tag(fields, "de:f:") {
        write!(
            out,
            "\tde:f:{:.6}",
            divergence(chunk.metrics.gap_compressed_nm, chunk.metrics.aln_len)
        )?;
    }
    if has_tag(fields, "df:i:") {
        write!(out, "\tdf:i:{}", chunk.metrics.nm)?;
    }
    write!(
        out,
        "\tzp:Z:{TOOL_NAME}\tzc:i:{idx}\tzl:i:{}\tzo:i:{}\tzm:Z:{}\tzs:i:{orig_q_start}\tze:i:{orig_q_end}\tzts:i:{orig_t_start}\tzte:i:{orig_t_end}\n",
        args.length,
        args.overlap,
        args.chunk_mode.as_str()
    )
}

fn is_alignment_derived_tag(tag: &str) -> bool {
    tag.starts_with("cg:Z:")
        || tag.starts_with("cs:Z:")
        || tag.starts_with("NM:i:")
        || tag.starts_with("dv:f:")
        || tag.starts_with("de:f:")
        || tag.starts_with("df:i:")
}

fn divergence(numer: u64, denom: u64) -> f64 {
    if denom == 0 {
        0.0
    } else {
        numer as f64 / denom as f64
    }
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
    let chunks = plan_chunks(fields, args, line_no)?;

    let mut chopped_bp = 0_u128;
    for (idx, chunk) in chunks.iter().enumerate() {
        emit_fragment(
            out, fields, args, idx as u64, q_start, q_end, t_start, t_end, chunk,
        )
        .map_err(|e| e.to_string())?;
        chopped_bp += (chunk.q_end - chunk.q_start) as u128;
    }
    Ok((chunks.len() as u64, chopped_bp))
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
        "comparison_id\ttool\tversion\tchop_length_bp\toverlap_bp\tchunk_mode\tkeep_tags\traw_records\tchopped_records\traw_query_bp\tchopped_query_bp"
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        args.comparison_id,
        TOOL_NAME,
        env!("CARGO_PKG_VERSION"),
        args.length,
        args.overlap,
        args.chunk_mode.as_str(),
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
        run_case_mode(input, length, keep_tags, ChunkMode::RowStart)
    }

    fn run_case_mode(
        input: &str,
        length: u64,
        keep_tags: bool,
        chunk_mode: ChunkMode,
    ) -> (String, Stats) {
        let args = Args {
            length,
            chunk_mode,
            keep_tags,
            comparison_id: "test".to_string(),
            ..Args::default()
        };
        let mut out = Vec::new();
        let stats = process(BufReader::new(input.as_bytes()), &mut out, &args).unwrap();
        (String::from_utf8(out).unwrap(), stats)
    }

    fn fields(line: &str) -> Vec<&str> {
        line.split('\t').collect()
    }

    #[test]
    fn splits_equal_match_chunks_and_recomputes_identity() {
        let input =
            "q\t100\t0\t25\t+\tt\t200\t100\t125\t25\t25\t60\tcg:Z:25=\tNM:i:0\tdv:f:0\tde:f:0\n";
        let (out, stats) = run_case(input, 10, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("q\t100\t0\t10\t+\tt\t200\t100\t110\t10\t10\t60\t"));
        assert!(lines[1].starts_with("q\t100\t10\t20\t+\tt\t200\t110\t120\t10\t10\t60\t"));
        assert!(lines[2].starts_with("q\t100\t20\t25\t+\tt\t200\t120\t125\t5\t5\t60\t"));
        assert!(lines[0].contains("cg:Z:10="));
        assert!(lines[0].contains("NM:i:0"));
        assert!(lines[0].contains("dv:f:0.000000"));
        assert!(lines[0].contains("de:f:0.000000"));
        assert!(lines[0].contains("zp:Z:pafchop-rs"));
        assert!(lines[0].contains("zm:Z:row-start"));
        assert_eq!(stats.raw_records, 1);
        assert_eq!(stats.chopped_records, 3);
        assert_eq!(stats.raw_query_bp, 25);
        assert_eq!(stats.chopped_query_bp, 25);
    }

    #[test]
    fn query_grid_uses_absolute_query_boundaries_for_shifted_rows() {
        let input = "\
q\t100\t7\t27\t+\tt1\t200\t100\t120\t20\t20\t60\tcg:Z:20=\n\
q\t100\t3\t23\t+\tt2\t200\t300\t320\t20\t20\t60\tcg:Z:20=\n";
        let (out, _) = run_case_mode(input, 10, false, ChunkMode::QueryGrid);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 6);
        let intervals: Vec<(&str, &str, &str)> = lines
            .iter()
            .map(|line| {
                let f = fields(line);
                (f[5], f[2], f[3])
            })
            .collect();
        assert_eq!(
            intervals,
            vec![
                ("t1", "7", "10"),
                ("t1", "10", "20"),
                ("t1", "20", "27"),
                ("t2", "3", "10"),
                ("t2", "10", "20"),
                ("t2", "20", "23"),
            ]
        );
        assert!(lines[1].contains("zm:Z:query-grid"));
        assert!(lines[4].starts_with("q\t100\t10\t20\t+\tt2\t200\t307\t317\t10\t10\t60\t"));
    }

    #[test]
    fn query_grid_reverse_strand_uses_absolute_boundaries() {
        let input = "q\t100\t7\t27\t-\tt\t200\t80\t100\t20\t20\t60\tcg:Z:20=\n";
        let (out, _) = run_case_mode(input, 10, false, ChunkMode::QueryGrid);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 3);
        let f0 = fields(lines[0]);
        assert_eq!(&f0[2..4], &["7", "10"]);
        assert_eq!(&f0[7..11], &["97", "100", "3", "3"]);
        let f1 = fields(lines[1]);
        assert_eq!(&f1[2..4], &["10", "20"]);
        assert_eq!(&f1[7..11], &["87", "97", "10", "10"]);
        let f2 = fields(lines[2]);
        assert_eq!(&f2[2..4], &["20", "27"]);
        assert_eq!(&f2[7..11], &["80", "87", "7", "7"]);
    }

    #[test]
    fn query_grid_clips_cigar_across_grid_boundary() {
        let input = "q\t100\t7\t27\t+\tt\t200\t50\t70\t14\t23\t60\tcg:Z:5=4X3I8=3D\tNM:i:10\tdv:f:0.1\tde:f:0.1\n";
        let (out, _) = run_case_mode(input, 10, false, ChunkMode::QueryGrid);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 3);

        let f0 = fields(lines[0]);
        assert_eq!(&f0[2..4], &["7", "10"]);
        assert_eq!(&f0[7..11], &["50", "53", "3", "3"]);
        assert!(lines[0].contains("cg:Z:3="));
        assert!(lines[0].contains("NM:i:0"));

        let f1 = fields(lines[1]);
        assert_eq!(&f1[2..4], &["10", "20"]);
        assert_eq!(&f1[7..11], &["53", "60", "3", "10"]);
        assert!(lines[1].contains("cg:Z:2=4X3I1="));
        assert!(lines[1].contains("NM:i:7"));

        let f2 = fields(lines[2]);
        assert_eq!(&f2[2..4], &["20", "27"]);
        assert_eq!(&f2[7..11], &["60", "70", "7", "10"]);
        assert!(lines[2].contains("cg:Z:7=3D"));
        assert!(lines[2].contains("NM:i:3"));
    }

    #[test]
    fn query_grid_output_matches_per_record_parallel_concatenation() {
        let input = "\
q1\t100\t7\t27\t+\tt1\t200\t50\t70\t20\t20\t60\tcg:Z:20=\n\
q2\t100\t3\t18\t-\tt2\t200\t80\t95\t15\t15\t60\tcg:Z:15=\n\
q3\t100\t11\t29\t+\tt3\t200\t100\t117\t12\t19\t60\tcg:Z:4=2X3I5=2D4=\tNM:i:7\n";
        let (sequential, sequential_stats) = run_case_mode(input, 10, false, ChunkMode::QueryGrid);

        let mut concatenated = String::new();
        let mut combined_stats = Stats::default();
        for line in input.lines() {
            let (chunked, stats) =
                run_case_mode(&format!("{line}\n"), 10, false, ChunkMode::QueryGrid);
            concatenated.push_str(&chunked);
            combined_stats.raw_records += stats.raw_records;
            combined_stats.chopped_records += stats.chopped_records;
            combined_stats.raw_query_bp += stats.raw_query_bp;
            combined_stats.chopped_query_bp += stats.chopped_query_bp;
        }

        assert_eq!(sequential, concatenated);
        assert_eq!(sequential_stats.raw_records, combined_stats.raw_records);
        assert_eq!(
            sequential_stats.chopped_records,
            combined_stats.chopped_records
        );
        assert_eq!(sequential_stats.raw_query_bp, combined_stats.raw_query_bp);
        assert_eq!(
            sequential_stats.chopped_query_bp,
            combined_stats.chopped_query_bp
        );
    }

    #[test]
    fn chunks_crossing_m_x_i_d_recompute_col10_col11() {
        let input = "q\t100\t0\t18\t+\tt\t200\t100\t117\t11\t20\t60\tcg:Z:4M2X3I5=2D2X2=\tcs:Z::4*ag*ct+aaa:5-tt*ga*tc:2\tNM:i:9\tdv:f:0.5\tde:f:0.5\tst:Z:keep\n";
        let (out, _) = run_case(input, 10, true);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 2);

        let f0 = fields(lines[0]);
        assert_eq!(&f0[2..4], &["0", "10"]);
        assert_eq!(&f0[7..11], &["100", "107", "5", "10"]);
        assert!(lines[0].contains("st:Z:keep"));
        assert!(lines[0].contains("cg:Z:4M2X3I1="));
        assert!(lines[0].contains("NM:i:5"));
        assert!(lines[0].contains("dv:f:0.500000"));
        assert!(lines[0].contains("de:f:0.300000"));

        let f1 = fields(lines[1]);
        assert_eq!(&f1[2..4], &["10", "18"]);
        assert_eq!(&f1[7..11], &["107", "117", "6", "10"]);
        assert!(lines[1].contains("cg:Z:4=2D2X2="));
        assert!(lines[1].contains("NM:i:4"));
        assert!(lines[1].contains("dv:f:0.400000"));
        assert!(lines[1].contains("de:f:0.300000"));
    }

    #[test]
    fn reverse_strand_uses_forward_target_coordinates() {
        let input = "q\t100\t0\t12\t-\tt\t200\t80\t95\t10\t15\t60\tcg:Z:5=3D5X2=\tNM:i:8\n";
        let (out, _) = run_case(input, 6, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 2);
        let f0 = fields(lines[0]);
        assert_eq!(&f0[2..4], &["0", "6"]);
        assert_eq!(&f0[7..11], &["86", "95", "5", "9"]);
        assert!(lines[0].contains("cg:Z:5=3D1X"));
        let f1 = fields(lines[1]);
        assert_eq!(&f1[2..4], &["6", "12"]);
        assert_eq!(&f1[7..11], &["80", "86", "2", "6"]);
        assert!(lines[1].contains("cg:Z:4X2="));
    }

    #[test]
    fn chunk_can_end_inside_cigar_operation() {
        let input = "q\t100\t0\t13\t+\tt\t200\t50\t63\t13\t13\t60\tcg:Z:13=\n";
        let (out, _) = run_case(input, 5, false);
        let cigars: Vec<&str> = out
            .lines()
            .map(|line| line.split('\t').find(|f| f.starts_with("cg:Z:")).unwrap())
            .collect();
        assert_eq!(cigars, vec!["cg:Z:5=", "cg:Z:5=", "cg:Z:3="]);
    }

    #[test]
    fn clips_cs_string_when_present() {
        let input = "q\t100\t0\t8\t+\tt\t200\t10\t17\t5\t9\t60\tcg:Z:3=1X2I2=1D\tcs:Z::3*ag+tt:2-c\tNM:i:4\n";
        let (out, _) = run_case(input, 5, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].contains("cg:Z:3=1X1I"));
        assert!(lines[0].contains("cs:Z::3*ag+t"));
        assert!(lines[1].contains("cg:Z:1I2=1D"));
        assert!(lines[1].contains("cs:Z:+t:2-c"));
    }

    #[test]
    fn refuses_to_chop_without_cigar_instead_of_interpolating_identity() {
        let input = "q\t100\t0\t25\t+\tt\t200\t100\t125\t20\t25\t60\n";
        let args = Args::default();
        let mut out = Vec::new();
        let err = process(BufReader::new(input.as_bytes()), &mut out, &args).unwrap_err();
        assert!(err.contains("cannot chop exactly without cg:Z"));
    }

    #[test]
    fn refuses_ambiguous_m_cigar_without_cs() {
        let input = "q\t100\t0\t8\t+\tt\t200\t100\t108\t8\t8\t60\tcg:Z:8M\n";
        let args = Args::default();
        let mut out = Vec::new();
        let err = process(BufReader::new(input.as_bytes()), &mut out, &args).unwrap_err();
        assert!(err.contains("ambiguous M operations"));
    }

    #[test]
    fn can_preserve_non_alignment_tags_when_requested() {
        let input =
            "q\t100\t0\t8\t+\tt\t200\t100\t108\t8\t8\t60\tcg:Z:8=\ttp:A:P\tst:Z:kept\tdv:f:0.0\n";
        let (out, stats) = run_case(input, 10, true);
        assert_eq!(stats.chopped_records, 1);
        assert!(out.contains("tp:A:P"));
        assert!(out.contains("st:Z:kept"));
        assert!(out.contains("cg:Z:8="));
        assert!(out.contains("dv:f:0.000000"));
    }

    #[test]
    fn supports_overlap() {
        let args = Args {
            length: 10,
            overlap: 2,
            comparison_id: "test".to_string(),
            ..Args::default()
        };
        let input = "q\t100\t0\t25\t+\tt\t200\t100\t125\t25\t25\t60\tcg:Z:25=\n";
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

    #[test]
    fn property_style_chunk_metrics_sum_to_original_without_overlap() {
        let cigars = ["7=2X3I5=2D4X6=", "2D10=1X1I1D7=", "3=2I2X2D3=4="];
        for cigar in cigars {
            let ops = parse_cigar(cigar, 1).unwrap();
            let q_len: u64 = ops
                .iter()
                .filter(|op| q_consumes(op.kind))
                .map(|op| op.len)
                .sum();
            let t_len: u64 = ops
                .iter()
                .filter(|op| t_consumes(op.kind))
                .map(|op| op.len)
                .sum();
            let mut original = ChunkMetrics::default();
            for op in &ops {
                add_metrics(&mut original, op.kind, op.len);
            }
            let input = format!(
                "q\t100\t0\t{q_len}\t+\tt\t200\t50\t{}\t{}\t{}\t60\tcg:Z:{cigar}\tNM:i:{}\n",
                50 + t_len,
                original.matches,
                original.aln_len,
                original.nm
            );
            let (out, _) = run_case(&input, 4, false);
            let mut sum_matches = 0_u64;
            let mut sum_aln = 0_u64;
            let mut sum_nm = 0_u64;
            for line in out.lines() {
                let f = fields(line);
                sum_matches += f[9].parse::<u64>().unwrap();
                sum_aln += f[10].parse::<u64>().unwrap();
                let nm = f.iter().find_map(|x| x.strip_prefix("NM:i:")).unwrap();
                sum_nm += nm.parse::<u64>().unwrap();
            }
            assert_eq!(sum_matches, original.matches, "{cigar}");
            assert_eq!(sum_aln, original.aln_len, "{cigar}");
            assert_eq!(sum_nm, original.nm, "{cigar}");
        }
    }
}
