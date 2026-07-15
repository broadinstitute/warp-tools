use anyhow::{Context, Result};
use rust_htslib::bcf::{
    header::{Header, HeaderRecord}, record::GenotypeAllele, Format, Read, IndexedReader, Writer
};
use std::collections::{HashMap, VecDeque};
use std::env;
use std::fs;

struct InputRec {
    rank: usize,
    pos: i64,
    end: i64,
    is_ref_block: bool,
    min_reps: Vec<(i64, Vec<u8>, Vec<u8>)>,
    gqs: Vec<i32>,
    pls: Vec<i32>,
    pl_stride: usize,
    gts: Vec<Option<usize>>,
}

fn get_minimal_representation<'a>(mut pos: i64, mut r: &'a [u8], mut a: &'a [u8]) -> (i64, &'a [u8], &'a [u8]) {
    if !a.is_empty() && a[0] == b'<' {
        return (pos, r, a);
    }
    while !r.is_empty() && !a.is_empty() && r.last() == a.last() {
        r = &r[..r.len() - 1];
        a = &a[..a.len() - 1];
    }
    while !r.is_empty() && !a.is_empty() && r.first() == a.first() {
        r = &r[1..];
        a = &a[1..];
        pos += 1;
    }
    (pos, r, a)
}

fn main() -> Result<()> {
    let mut args = env::args().skip(1);
    if args.len() < 4 {
        eprintln!("Usage: extract_bubble_PLs <gvcf|joint> <panel.bcf> <input.vcf.gz> <output.bcf> [--region chr:start-end] [--samples list.txt] [--window 15000] [--cap-pl 30] [--scale-pl 5.0] [--threads 4]");
        std::process::exit(1);
    }

    let mode = args.next().unwrap();
    let is_joint = match mode.as_str() {
        "gvcf" => false,
        "joint" => true,
        _ => {
            eprintln!("Error: First argument must be strictly 'gvcf' or 'joint'.");
            std::process::exit(1);
        }
    };

    // --- MINIMAL CHANGE: Silent Custom Index Adapter ---
    let handle_custom_index = |arg: String| -> Result<String> {
        if let Some(pos) = arg.find("##idx##") {
            let vcf_path = &arg[..pos];
            let idx_path = &arg[pos + 7..];

            let expected_idx = if idx_path.ends_with(".tbi") {
                format!("{}.tbi", vcf_path)
            } else {
                format!("{}.csi", vcf_path)
            };

            if !std::path::Path::new(&expected_idx).exists() {
                std::os::unix::fs::symlink(idx_path, &expected_idx)
                    .with_context(|| format!("Failed to symlink custom index {} to {}", idx_path, expected_idx))?;
            }
            Ok(vcf_path.to_string())
        } else {
            Ok(arg)
        }
    };

    let panel_path = handle_custom_index(args.next().unwrap())?;
    let input_path = handle_custom_index(args.next().unwrap())?;
    let out_path = args.next().unwrap();
    // --------------------------------------------------

    let mut region_str: Option<String> = None;
    let mut samples_file: Option<String> = None;
    let mut window_size = 15000i64;
    let mut cap_pl = 30i32;
    let mut scale_pl = 5.0f64;
    let mut threads = 4usize;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--region" => region_str = Some(args.next().unwrap()),
            "--samples" => samples_file = Some(args.next().unwrap()),
            "--window" => window_size = args.next().unwrap().parse().unwrap(),
            "--cap-pl" => cap_pl = args.next().unwrap().parse().unwrap(),
            "--scale-pl" => scale_pl = args.next().unwrap().parse().unwrap(),
            "--threads" => threads = args.next().unwrap().parse().unwrap(),
            _ => {}
        }
    }

    eprintln!("[DEBUG] Opening panel BCF using IndexedReader for sharding...");
    let mut panel_reader = IndexedReader::from_path(panel_path).context("Failed to open Panel BCF. Is it indexed?")?;
    panel_reader.set_threads(threads).context("Failed to set threads on panel reader")?;
    let panel_header = panel_reader.header().clone();

    // REGION PARSING LOGIC
    let mut target_rid = None;
    let mut target_start = 0i64;
    let mut target_end = i64::MAX;

    if let Some(reg) = &region_str {
        let parts: Vec<&str> = reg.split(':').collect();
        let chrom = parts[0];
        target_rid = Some(panel_header.name2rid(chrom.as_bytes()).context("Region contig not found in panel dictionary")? as usize);

        if parts.len() > 1 {
            let pos_parts: Vec<&str> = parts[1].split('-').collect();
            target_start = pos_parts[0].parse::<i64>().context("Invalid region start")? - 1;
            if pos_parts.len() > 1 && !pos_parts[1].trim().is_empty() {
                target_end = pos_parts[1].parse::<i64>().context("Invalid region end")? - 1;
            } else {
                target_end = target_start;
            }
        }

        let end_opt = if target_end == i64::MAX { None } else { Some(target_end as u64) };
        panel_reader.fetch(target_rid.unwrap() as u32, target_start as u64, end_opt).context("Failed to fetch region on panel")?;
        eprintln!("[DEBUG] Region target established. Contig: {}, bounds: [{}, {}]", chrom, target_start, target_end);
    }

    eprintln!("[DEBUG] Opening input VCF/gVCF using IndexedReader...");
    let mut input_reader = IndexedReader::from_path(input_path).context("Failed to open input VCF. Is it indexed?")?;
    input_reader.set_threads(threads).context("Failed to set threads on input reader")?;
    let input_header = input_reader.header().clone();

    // SAMPLE PARSING LOGIC
    let input_samples: Vec<&[u8]> = input_header.samples();
    let mut sample_mapping: Vec<usize> = Vec::new();

    if let Some(path) = &samples_file {
        eprintln!("[DEBUG] Parsing requested samples from {}...", path);
        let content = fs::read_to_string(path).context("Failed to read samples file")?;
        for line in content.lines() {
            let s = line.trim();
            if s.is_empty() { continue; }
            if let Some(idx) = input_samples.iter().position(|&in_s| in_s == s.as_bytes()) {
                sample_mapping.push(idx);
            } else {
                eprintln!("[WARN] Requested sample '{}' not found in input VCF. Skipping.", s);
            }
        }
    } else {
        sample_mapping = (0..input_samples.len()).collect();
    }

    let out_n_samples = sample_mapping.len();
    if out_n_samples == 0 {
        eprintln!("[ERROR] No overlapping samples found to process.");
        std::process::exit(1);
    }

    eprintln!("[DEBUG] Building unified dictionary from input header...");
    let mut chrom_rank = HashMap::new();
    for record in input_header.header_records() {
        if let HeaderRecord::Contig { values, .. } = record {
            if let Some(id) = values.get("ID") {
                if let Ok(rid) = input_header.name2rid(id.as_bytes()) {
                    chrom_rank.insert(id.as_bytes().to_vec(), rid as usize);
                }
            }
        }
    }

    eprintln!("[DEBUG] Subsetting template and pushing FORMAT tags...");
    let empty_samples: &[&[u8]] = &[];
    let mut out_header = Header::from_template_subset(&panel_header, empty_samples).unwrap();
    out_header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    out_header.push_record(b"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods\">");

    for &idx in &sample_mapping {
        out_header.push_sample(input_samples[idx]);
    }

    eprintln!("[DEBUG] Creating writer...");
    let mut out_writer = Writer::from_path(out_path, &out_header, false, Format::Bcf).unwrap();
    out_writer.set_threads(threads).context("Failed to set threads on output writer")?;
    let mut out_record = out_writer.empty_record();
    let mut g_rec = input_reader.empty_record();

    let mut input_buffer: VecDeque<InputRec> = VecDeque::new();
    let mut input_eof = false;
    let mut records_processed = 0usize;

    let process_pl = |v: i32| -> i32 {
        if v == i32::MIN { return v; }
        let scaled = ((v as f64) / scale_pl) as i32;
        if scaled > cap_pl { cap_pl } else { scaled }
    };

    eprintln!("[INFO] Mode={}, Output Samples={}, Window={}, Cap={}, Scale={}, Threads={}",
        mode, out_n_samples, window_size, cap_pl, scale_pl, threads);

    let mut current_panel_rid = usize::MAX;
    let mut current_p_rank = usize::MAX;
    let mut current_panel_chrom_str = String::new();

    let mut out_gt = vec![GenotypeAllele::UnphasedMissing; out_n_samples * 2];
    let mut out_pl = vec![0i32; out_n_samples * 3];

    for p_rec_res in panel_reader.records() {
        let p_rec = p_rec_res?;
        let p_pos = p_rec.pos();

        if target_rid.is_some() && (p_pos < target_start || p_pos > target_end) {
            continue;
        }

        records_processed += 1;
        out_record.clear();

        let p_rid = p_rec.rid().unwrap() as usize;

        if p_rid != current_panel_rid {
            current_panel_rid = p_rid;
            let p_chrom_bytes = panel_header.rid2name(p_rid as u32).unwrap();
            current_p_rank = *chrom_rank.get(p_chrom_bytes).unwrap_or(&usize::MAX);

            current_panel_chrom_str = String::from_utf8_lossy(p_chrom_bytes).into_owned();
            eprintln!("\n[DEBUG] Panel switched to Contig: {} (Mapped to Input Rank: {})", current_panel_chrom_str, current_p_rank);

            input_buffer.clear();
            input_eof = false;

            if current_p_rank != usize::MAX {
                // DYNAMIC SHARD SEEK: Only apply the massive 100kb backward padding for gVCFs.
                // Joint VCFs just need the standard window size.
                let padding = if is_joint { 0 } else { 100_000 };
                let fetch_start = (p_pos - window_size - padding).max(0) as u64;

                eprintln!("[DEBUG] Synchronizing input reader to {}:{}...", current_panel_chrom_str, fetch_start);

                if let Err(e) = input_reader.fetch(current_p_rank as u32, fetch_start, None) {
                    eprintln!("[DEBUG-WARN] Could not fetch contig: {}", e);
                    input_eof = true;
                }
            }
        }

        let p_alleles = p_rec.alleles();
        if p_alleles.len() < 2 { continue; }

        if records_processed == 1 {
            eprintln!("[Progress] Buffer initialized! Processing first record at {}:{}", current_panel_chrom_str, p_pos + 1);
        } else if records_processed % 100_000 == 0 {
            eprintln!("[Progress] Scanned {} panel sites... currently at {}:{}", records_processed, current_panel_chrom_str, p_pos + 1);
        }

        let min_p = get_minimal_representation(p_pos, p_alleles[0], p_alleles[1]);
        let mut matched = false;

        out_gt.fill(GenotypeAllele::UnphasedMissing);
        out_pl.fill(0);

        // --- 1. ADVANCE BUFFER ---
        while !input_eof {
            if let Some(back) = input_buffer.back() {
                if back.rank > current_p_rank || (back.rank == current_p_rank && back.pos > p_pos + window_size + 50) {
                    break;
                }
            }

            match input_reader.read(&mut g_rec) {
                Some(Ok(_)) => {
                    let g_rid = g_rec.rid().unwrap();
                    let g_chrom_bytes = input_header.rid2name(g_rid).unwrap();
                    let g_rank = *chrom_rank.get(g_chrom_bytes).unwrap_or(&usize::MAX);

                    if g_rank == usize::MAX { continue; }

                    let pos = g_rec.pos();
                    let alleles = g_rec.alleles();
                    let is_ref_block = alleles.len() > 1 && alleles[1] == b"<NON_REF>";

                    let end = if is_ref_block {
                        if let Ok(Some(end_vals)) = g_rec.info(b"END").integer() { end_vals[0] as i64 - 1 } else { pos }
                    } else { pos + alleles[0].len() as i64 - 1 };

                    if end < p_pos - window_size - 50_000 { continue; }

                    let mut min_reps = Vec::new();
                    let mut gts = vec![None; out_n_samples * 2];
                    let mut gqs = vec![0; out_n_samples];

                    // FIXED: Mathematically strict PL stride calculation
                    let n_alleles = alleles.len();
                    let pl_stride = (n_alleles * (n_alleles + 1)) / 2;
                    let mut pls = vec![i32::MIN; out_n_samples * pl_stride];

                    if is_ref_block {
                        if let Ok(gq_format) = g_rec.format(b"GQ").integer() {
                            for (out_idx, &in_idx) in sample_mapping.iter().enumerate() {
                                gqs[out_idx] = match gq_format.get(in_idx) {
                                    Some(slice) if !slice.is_empty() => slice[0],
                                    _ => 0,
                                };
                            }
                        }
                    } else {
                        for alt in &alleles[1..] {
                            if *alt == b"<NON_REF>" {
                                min_reps.push((0, vec![], vec![]));
                            } else {
                                let min_g = get_minimal_representation(pos, alleles[0], *alt);
                                min_reps.push((min_g.0, min_g.1.to_vec(), min_g.2.to_vec()));
                            }
                        }

                        // Safely extract LPL or PL, strictly bound to the computed pl_stride
                        let mut pl_found = false;
                        if let Ok(lpl_fmt) = g_rec.format(b"LPL").integer() {
                            for (out_idx, &in_idx) in sample_mapping.iter().enumerate() {
                                if let Some(slice) = lpl_fmt.get(in_idx) {
                                    let base = out_idx * pl_stride;
                                    for i in 0..slice.len().min(pl_stride) {
                                        pls[base + i] = slice[i];
                                    }
                                }
                            }
                            pl_found = true;
                        }

                        if !pl_found {
                            if let Ok(pl_fmt) = g_rec.format(b"PL").integer() {
                                for (out_idx, &in_idx) in sample_mapping.iter().enumerate() {
                                    if let Some(slice) = pl_fmt.get(in_idx) {
                                        let base = out_idx * pl_stride;
                                        for i in 0..slice.len().min(pl_stride) {
                                            pls[base + i] = slice[i];
                                        }
                                    }
                                }
                            }
                        }

                        if let Ok(gt_fmt) = g_rec.genotypes() {
                            for (out_idx, &in_idx) in sample_mapping.iter().enumerate() {
                                let gt_slice = gt_fmt.get(in_idx);
                                if gt_slice.len() == 2 {
                                    gts[out_idx * 2] = gt_slice[0].index().map(|idx| idx as usize);
                                    gts[out_idx * 2 + 1] = gt_slice[1].index().map(|idx| idx as usize);
                                }
                            }
                        }
                    }

                    input_buffer.push_back(InputRec {
                        rank: g_rank, pos, end, is_ref_block,
                        min_reps, gqs, pls, pl_stride, gts
                    });
                },
                Some(Err(_)) | None => { input_eof = true; }
            }
        }

        // --- 2. PRUNE BUFFER ---
        while let Some(front) = input_buffer.front() {
            if front.rank < current_p_rank || (front.rank == current_p_rank && front.end < p_pos - window_size - 50) {
                input_buffer.pop_front();
            } else { break; }
        }

        // --- 3. SEARCH BUFFER (Variant Match) ---
        for g in input_buffer.iter() {
            if g.rank != current_p_rank || g.is_ref_block { continue; }

            for (k_idx, min_g) in g.min_reps.iter().enumerate() {
                if min_g.1.is_empty() && min_g.2.is_empty() { continue; } // skip only the <NON_REF> placeholder, not real insertions
                let k = k_idx + 1;

                if min_g.0 == min_p.0 && min_g.1 == min_p.1 && min_g.2 == min_p.2 {
                    let idx0 = 0;
                    let idx1 = k * (k + 1) / 2;
                    let idx2 = idx1 + k;

                    for s in 0..out_n_samples {
                        let map_gt = |idx: Option<usize>| match idx {
                            Some(0) => GenotypeAllele::Unphased(0),
                            Some(v) if v == k => GenotypeAllele::Unphased(1),
                            _ => GenotypeAllele::UnphasedMissing,
                        };
                        out_gt[s * 2] = map_gt(g.gts[s * 2]);
                        out_gt[s * 2 + 1] = map_gt(g.gts[s * 2 + 1]);

                        if g.pl_stride > idx2 {
                            let base = s * g.pl_stride;
                            out_pl[s * 3] = process_pl(g.pls[base + idx0]);
                            out_pl[s * 3 + 1] = process_pl(g.pls[base + idx1]);
                            out_pl[s * 3 + 2] = process_pl(g.pls[base + idx2]);
                        }
                    }
                    matched = true;
                    break;
                }
            }
            if matched { break; }
        }

        // --- 4. FALLBACK LOGIC ---
        if !matched {
            if is_joint {
                // DYNAMIC FILTER: If the bubble isn't in the joint VCF, completely
                // drop the row to save disk space and GLIMPSE2 compute time.
                continue;
            } else {
                for g in input_buffer.iter() {
                    if g.rank != current_p_rank || !g.is_ref_block { continue; }
                    if g.pos <= min_p.0 && g.end >= min_p.0 + min_p.1.len() as i64 - 1 {
                        for s in 0..out_n_samples {
                            let gq = if g.gqs[s] == i32::MIN { 0 } else { g.gqs[s] };
                            out_gt[s * 2] = GenotypeAllele::Unphased(0);
                            out_gt[s * 2 + 1] = GenotypeAllele::Unphased(0);
                            out_pl[s * 3] = 0;
                            out_pl[s * 3 + 1] = process_pl(gq);
                            out_pl[s * 3 + 2] = process_pl(gq * 2);
                        }
                        break;
                    }
                }
            }
        }

        out_record.set_rid(Some(p_rid as u32));
        out_record.set_pos(p_pos);
        out_record.set_alleles(&[p_alleles[0], p_alleles[1]])?;
        if let Ok(Some(info_ids)) = p_rec.info(b"ID").string() {
            out_record.set_id(info_ids[0])?;
        } else {
            out_record.set_id(b".")?;
        }
        out_record.push_genotypes(&out_gt)?;
        out_record.push_format_integer(b"PL", &out_pl)?;
        out_writer.write(&out_record)?;
    }

    eprintln!("[INFO] Loop complete. Flushing background compression threads to disk...");
    drop(out_writer);
    println!("[INFO] Preprocessed imputation BCF generated successfully.");
    std::process::exit(0);
}
