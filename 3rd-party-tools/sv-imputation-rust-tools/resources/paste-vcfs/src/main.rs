use clap::Parser;
use rust_htslib::bcf::{self, header::{Header, TagType}, Read};
use rust_htslib::htslib;
use std::ffi::CString;
use std::os::raw::{c_int, c_void};
use std::ptr;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about = "Pastes BCF/VCF.GZ across identical sites; index-seek for regions")]
struct Args {
    #[arg(short, long)]
    region: Option<String>,

    #[arg(short, long, default_value_t = 1)]
    threads: usize,

    #[arg(long, value_delimiter = ',')]
    info: Vec<String>,

    #[arg(long, value_delimiter = ',')]
    format: Vec<String>,

    #[arg(short, long)]
    output: String,

    #[arg(required = true)]
    inputs: Vec<String>,
}

struct TagSpec {
    bytes: Vec<u8>,
    ty: TagType,
}

/// Lightweight per-file region reader: index-seek + one-record-at-a-time
/// streaming via hts_itr_next. No synced-reader buffering.
struct RegionReader {
    fp: *mut htslib::htsFile,
    hdr: *mut htslib::bcf_hdr_t,
    idx: *mut htslib::hts_idx_t,
    itr: *mut htslib::hts_itr_t,
}

impl RegionReader {
    // `tid` is the contig's numeric id, resolved once by the caller via the
    // safe name2rid API. This avoids the bcf_hdr_name2id macro, which bindgen
    // does not emit (only the underlying bcf_hdr_id2int is exported). Inputs
    // share contig order, so one tid serves every file.
    fn new(path: &str, tid: i32, beg: i64, end: i64) -> Self {
        let c_path = CString::new(path).expect("path has interior NUL");
        let mode = CString::new("r").unwrap();
        unsafe {
            let fp = htslib::hts_open(c_path.as_ptr(), mode.as_ptr());
            assert!(!fp.is_null(), "failed to open {}", path);

            let hdr = htslib::bcf_hdr_read(fp);
            assert!(!hdr.is_null(), "failed to read header from {}", path);

            // bcf_index_load(fn) macro == hts_idx_load(fn, HTS_FMT_CSI)
            let idx = htslib::hts_idx_load(c_path.as_ptr(), htslib::HTS_FMT_CSI as c_int);
            assert!(!idx.is_null(), "no .csi index for {} (run `bcftools index`)", path);

            // bcf_itr_queryi(idx, tid, beg, end) macro
            //   == hts_itr_query(idx, tid, beg, end, bcf_readrec)
            let itr = htslib::hts_itr_query(idx, tid, beg, end, Some(htslib::bcf_readrec));
            assert!(!itr.is_null(), "failed to query region (tid {}) in {}", tid, path);

            RegionReader { fp, hdr, idx, itr }
        }
    }

    /// Reads the next record in the region into `record`.
    /// Returns Some(()) on success, None at end of region.
    fn read_into(&mut self, record: &mut bcf::Record) -> Option<()> {
        unsafe {
            // bcf_itr_next(fp, itr, r) macro
            //   == hts_itr_next(fp->fp.bgzf, itr, r, 0)
            let bgzf = (*self.fp).fp.bgzf;
            let ret = htslib::hts_itr_next(
                bgzf,
                self.itr,
                record.inner() as *const _ as *mut c_void, // inner() accessor + const-cast for FFI write
                ptr::null_mut(),
            );
            if ret >= 0 { Some(()) } else { None }
        }
    }
}

impl Drop for RegionReader {
    fn drop(&mut self) {
        unsafe {
            if !self.itr.is_null() { htslib::hts_itr_destroy(self.itr); }
            if !self.idx.is_null() { htslib::hts_idx_destroy(self.idx); }
            if !self.hdr.is_null() { htslib::bcf_hdr_destroy(self.hdr); }
            if !self.fp.is_null() { htslib::hts_close(self.fp); }
        }
    }
}

/// Unifies whole-file streaming and region-iterator reading.
enum Source {
    Whole(bcf::Reader),
    Region(RegionReader),
}

impl Source {
    fn read_into(&mut self, rec: &mut bcf::Record) -> Option<()> {
        match self {
            Source::Whole(r) => match r.read(rec) {
                Some(Ok(())) => Some(()),
                Some(Err(e)) => panic!("Error reading record: {}", e),
                None => None,
            },
            Source::Region(r) => r.read_into(rec),
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Setup readers: used only to build the merged header and resolve tag
    // types. They stay open to provide a header for rid2name in progress.
    let setup_readers: Vec<bcf::Reader> = args.inputs.iter()
        .map(|path| bcf::Reader::from_path(path).expect("Failed to open input file"))
        .collect();

    let num_files = setup_readers.len();
    if num_files == 0 {
        panic!("No input files provided.");
    }

    let mut merged_header = Header::from_template(setup_readers[0].header());
    for reader in setup_readers.iter().skip(1) {
        for sample in reader.header().samples() {
            merged_header.push_sample(sample);
        }
    }

    let mut writer = bcf::Writer::from_path(&args.output, &merged_header, false, bcf::Format::Bcf)?;
    if args.threads > 1 {
        writer.set_threads(args.threads).expect("Failed to set writer threads");
    }

    let info_specs: Vec<TagSpec> = args.info.iter().map(|s| {
        let bytes = s.as_bytes().to_vec();
        let ty = setup_readers[0].header().info_type(&bytes).unwrap().0;
        TagSpec { bytes, ty }
    }).collect();

    let format_specs: Vec<TagSpec> = args.format.iter().map(|s| {
        let bytes = s.as_bytes().to_vec();
        let ty = if bytes == b"GT" {
            TagType::Integer
        } else {
            setup_readers[0].header().format_type(&bytes).unwrap().0
        };
        TagSpec { bytes, ty }
    }).collect();

    // Extract the regions to process (split by comma)
    let regions_to_process: Vec<Option<String>> = match &args.region {
        Some(r_str) => r_str.split(',').map(|s| Some(s.to_string())).collect(),
        None => vec![None],
    };

    // 1. Create records tied to their TRUE original headers, NOT the writer.
    // We do this BEFORE we consume setup_readers.
    let mut base_record = setup_readers[0].empty_record();
    let mut side_records: Vec<bcf::Record> = setup_readers.iter().skip(1).map(|r| r.empty_record()).collect();

    // new_record correctly belongs to the writer's merged header.
    let mut new_record = writer.empty_record();

    let mut int_buffer: Vec<i32> = Vec::new();
    let mut float_buffer: Vec<f32> = Vec::new();
    let mut string_storage: Vec<Vec<u8>> = Vec::new();

    let mut record_count: usize = 0;
    let start_time = Instant::now();

    // Wrap setup_readers in an Option so we can safely give ownership to Source::Whole later
    // or keep the readers alive indefinitely for Region loop
    let mut opt_setup_readers = Some(setup_readers);

    for reg_opt in &regions_to_process {
        // Parse region into iterator bounds + in-memory guard bounds.
        let mut req_rid: Option<u32> = None;
        let mut req_start = i64::MIN;
        let mut req_end = i64::MAX;
        let mut region_parts: Option<(String, i64, i64)> = None; // (chrom, iter_beg, iter_end)

        if let Some(r) = reg_opt {
            // Borrow the setup readers safely to query the header
            let readers_ref = opt_setup_readers.as_ref().expect("Readers consumed prematurely");
            let (chrom, coords) = r.split_once(':').unwrap_or((r, ""));
            let rid = readers_ref[0].header().name2rid(chrom.as_bytes())
                .unwrap_or_else(|_| panic!("Chromosome {} not found in header", chrom));
            req_rid = Some(rid);

            let mut start: u64 = 0;
            if !coords.is_empty() {
                let parts: Vec<&str> = coords.split('-').collect();
                if !parts.is_empty() && !parts[0].is_empty() {
                    let s = parts[0].replace(",", "").parse::<u64>().unwrap_or(1);
                    start = s.saturating_sub(1);
                    req_start = start as i64;
                }
                if parts.len() > 1 && !parts[1].is_empty() {
                    let e = parts[1].replace(",", "").parse::<u64>().unwrap_or(u64::MAX);
                    req_end = e.saturating_sub(1) as i64;
                } else if parts.len() == 1 {
                    req_end = start as i64;
                }
            }

            let iter_beg = if req_start == i64::MIN { 0 } else { req_start };
            let iter_end = if req_end == i64::MAX { i64::MAX } else { req_end + 1 }; // half-open
            region_parts = Some((chrom.to_string(), iter_beg, iter_end));
        }

        // 2. Build the read sources.
        let mut sources: Vec<Source> = Vec::with_capacity(num_files);

        if let Some((ref _chrom, beg, end)) = region_parts {
            let tid = req_rid.expect("region implies a resolved rid") as i32;
            for p in &args.inputs {
                sources.push(Source::Region(RegionReader::new(p, tid, beg, end)));
            }
        } else {
            // In Whole mode, Source::Whole safely takes ownership of the readers.
            let readers = opt_setup_readers.take().expect("Cannot read whole file multiple times");
            for r in readers {
                sources.push(Source::Whole(r));
            }
        }

        loop {
            match sources[0].read_into(&mut base_record) {
                Some(()) => {}
                None => break,
            }
            for i in 1..num_files {
                match sources[i].read_into(&mut side_records[i - 1]) {
                    Some(()) => {}
                    None => panic!("FATAL: Side file ran out of records prematurely!"),
                }
            }

            let base_rid = base_record.rid();
            let pos = base_record.pos();

            // In-memory guard. Redundant in region mode (the iterator already
            // restricts), harmless otherwise. Drops us out of inner loop when the region is exhausted.
            if let Some(want) = req_rid {
                if base_rid != Some(want) { continue; }
                if pos < req_start { continue; }
                if pos > req_end { break; }
            }

            let rid = base_rid.unwrap_or(0);

            record_count += 1;
            if record_count % 10_000 == 0 {
                // setup_readers may have been consumed in whole-file mode; use a
                // surviving header. In region mode setup_readers is still alive,
                // but to be uniform we look the name up via the writer header.
                let chrom_str = writer.header().rid2name(rid)
                    .map(|b| String::from_utf8_lossy(b).into_owned())
                    .unwrap_or_else(|_| "unknown".to_string());
                eprintln!(
                    "[Progress] Pasted {} records in {:.2?} | Current: {}:{}",
                    record_count, start_time.elapsed(), chrom_str, pos + 1
                );
            }

            let base_id = base_record.id();
            let base_alleles = base_record.alleles();
            for side_record in &side_records {
                if base_rid != side_record.rid() || pos != side_record.pos() {
                    panic!("FATAL: Position mismatch!");
                }
                if base_id != side_record.id() {
                    panic!("FATAL: ID mismatch!");
                }
                if base_alleles != side_record.alleles() {
                    panic!("FATAL: REF/ALT mismatch!");
                }
            }

            new_record.clear();
            new_record.set_rid(base_rid);
            new_record.set_pos(pos);
            new_record.set_id(&base_id)?;
            new_record.set_alleles(&base_alleles)?;
            new_record.set_qual(base_record.qual());

            for spec in &info_specs {
                let tag = spec.bytes.as_slice();
                match spec.ty {
                    TagType::Integer => {
                        if let Some(val) = base_record.info(tag).integer().unwrap_or(None).as_ref().map(|b| &**b) {
                            new_record.push_info_integer(tag, val)?;
                        }
                    }
                    TagType::Float => {
                        if let Some(val) = base_record.info(tag).float().unwrap_or(None).as_ref().map(|b| &**b) {
                            new_record.push_info_float(tag, val)?;
                        }
                    }
                    TagType::String => {
                        if let Some(val) = base_record.info(tag).string().unwrap_or(None).as_ref().map(|b| &**b) {
                            new_record.push_info_string(tag, val)?;
                        }
                    }
                    TagType::Flag => {
                        if base_record.info(tag).flag().unwrap_or(false) {
                            new_record.push_info_flag(tag)?;
                        }
                    }
                }
            }

            for spec in &format_specs {
                let tag = spec.bytes.as_slice();
                let tag_name = String::from_utf8_lossy(tag);

                // Closure to format the position, pinpoint the file, and panic
                let panic_mismatch = |missing_in_base: bool, side_idx: usize| -> ! {
                    let chrom_str = writer.header().rid2name(rid)
                        .map(|b| String::from_utf8_lossy(b).into_owned())
                        .unwrap_or_else(|_| "unknown".to_string());

                    let base_file = &args.inputs[0];
                    let side_file = &args.inputs[side_idx];

                    if missing_in_base {
                        panic!(
                            "FATAL: FORMAT tag '{}' mismatch at {}:{}. Tag is missing in base file '{}' but present in '{}'.",
                            tag_name, chrom_str, pos + 1, base_file, side_file
                        );
                    } else {
                        panic!(
                            "FATAL: FORMAT tag '{}' mismatch at {}:{}. Tag is present in base file '{}' but missing in '{}'.",
                            tag_name, chrom_str, pos + 1, base_file, side_file
                        );
                    }
                };

                match spec.ty {
                    TagType::Integer => {
                        int_buffer.clear();
                        if let Ok(vals) = base_record.format(tag).integer() {
                            for v in vals.iter() { int_buffer.extend_from_slice(v); }

                            for (i, side_record) in side_records.iter_mut().enumerate() {
                                if let Ok(side_vals) = side_record.format(tag).integer() {
                                    for v in side_vals.iter() { int_buffer.extend_from_slice(v); }
                                } else {
                                    panic_mismatch(false, i + 1);
                                }
                            }
                            new_record.push_format_integer(tag, &int_buffer)?;
                        } else {
                            // Base is missing the tag. Ensure sides are also missing it.
                            for (i, side_record) in side_records.iter_mut().enumerate() {
                                if side_record.format(tag).integer().is_ok() {
                                    panic_mismatch(true, i + 1);
                                }
                            }
                        }
                    }
                    TagType::Float => {
                        float_buffer.clear();
                        if let Ok(vals) = base_record.format(tag).float() {
                            for v in vals.iter() { float_buffer.extend_from_slice(v); }

                            for (i, side_record) in side_records.iter_mut().enumerate() {
                                if let Ok(side_vals) = side_record.format(tag).float() {
                                    for v in side_vals.iter() { float_buffer.extend_from_slice(v); }
                                } else {
                                    panic_mismatch(false, i + 1);
                                }
                            }
                            new_record.push_format_float(tag, &float_buffer)?;
                        } else {
                            for (i, side_record) in side_records.iter_mut().enumerate() {
                                if side_record.format(tag).float().is_ok() {
                                    panic_mismatch(true, i + 1);
                                }
                            }
                        }
                    }
                    TagType::String => {
                        let mut n = 0usize;
                        if let Ok(vals) = base_record.format(tag).string() {
                            for v in vals.iter() {
                                if n == string_storage.len() { string_storage.push(Vec::new()); }
                                string_storage[n].clear();
                                string_storage[n].extend_from_slice(v);
                                n += 1;
                            }

                            for (i, side_record) in side_records.iter_mut().enumerate() {
                                if let Ok(side_vals) = side_record.format(tag).string() {
                                    for v in side_vals.iter() {
                                        if n == string_storage.len() { string_storage.push(Vec::new()); }
                                        string_storage[n].clear();
                                        string_storage[n].extend_from_slice(v);
                                        n += 1;
                                    }
                                } else {
                                    panic_mismatch(false, i + 1);
                                }
                            }
                            let slices: Vec<&[u8]> = string_storage[..n].iter().map(|v| v.as_slice()).collect();
                            new_record.push_format_string(tag, &slices)?;
                        } else {
                            for (i, side_record) in side_records.iter_mut().enumerate() {
                                if side_record.format(tag).string().is_ok() {
                                    panic_mismatch(true, i + 1);
                                }
                            }
                        }
                    }
                    _ => panic!("Unsupported FORMAT type for tag {}", tag_name),
                }
            }

            writer.write(&new_record)?;
        }
    }

    Ok(())
}
