use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, BufWriter};
use std::time::Instant;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

struct Record {
    raf: Option<String>,
    af: Option<String>,
    info_val: Option<String>,
    atomic_ids: Vec<String>,
}

#[derive(Clone)]
struct PhasedDist {
    p0: f32,
    p1: f32,
}

fn smart_open(filename: &str) -> Box<dyn BufRead> {
    let file = File::open(filename).expect("Cannot open file");
    if filename.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

fn format_float(val: f32) -> String {
    let s = format!("{:.3}", val);
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

fn process_group(
    group_lines: &[(String, String)],
    id_buffer: &HashMap<String, (u32, String, String, String)>,
    max_alleles: usize,
    out_handle: &mut impl Write,
) {
    if group_lines.is_empty() { return; }

    let first_line_fields: Vec<&str> = group_lines[0].0.trim_end().split('\t').collect();

    if first_line_fields.len() <= 9 {
        panic!("Error: VCF does not contain sample columns. This script requires sample Genotype (GT) and Probability (GP) columns to project joint distributions.");
    }

    let num_samples = first_line_fields.len() - 9;
    let chrom = first_line_fields[0].to_string();
    let num_alleles = group_lines.len();

    let mut records: Vec<Record> = Vec::with_capacity(num_alleles);
    let mut all_atomic_ids = HashSet::new();

    let mut hap_probs: Vec<Vec<(f32, f32)>> = vec![vec![(0.0, 0.0); num_alleles]; num_samples];

    for (a, (line, site_info)) in group_lines.iter().enumerate() {
        let fields: Vec<&str> = line.trim_end().split('\t').collect();

        let mut raf = None;
        let mut af = None;
        let mut info_val = None;

        for item in fields[7].split(';') {
            if let Some(v) = item.strip_prefix("RAF=") { raf = Some(v.to_string()); }
            else if let Some(v) = item.strip_prefix("AF=") { af = Some(v.to_string()); }
            else if let Some(v) = item.strip_prefix("INFO=") { info_val = Some(v.to_string()); }
        }

        let mut atomic_ids = Vec::new();
        for item in site_info.split(';') {
            if let Some(id_str) = item.strip_prefix("ID=") {
                let replaced_id = id_str.replace(',', ":");
                for j in replaced_id.split(':').map(|s| s.trim()) {
                    if id_buffer.contains_key(j) {
                        atomic_ids.push(j.to_string());
                        all_atomic_ids.insert(j.to_string());
                    } else {
                        panic!("Error: Variant ID '{}' not found in the ID buffer. Ensure your biallelic VCF contains this ID and the window size is large enough.", j);
                    }
                }
                break;
            }
        }

        records.push(Record {
            raf, af, info_val, atomic_ids,
        });

        let fmt: Vec<&str> = fields[8].split(':').collect();
        let gt_idx = fmt.iter().position(|&x| x == "GT");
        let gp_idx = fmt.iter().position(|&x| x == "GP");

        for s in 0..num_samples {
            let sample_data = fields[9 + s];
            let mut gt_val = "0|0";
            let mut gp1 = 0.0;
            let mut gp2 = 0.0;

            if sample_data != "." {
                let mut split_iter = sample_data.split(':');
                let mut current_idx = 0;

                while let Some(val) = split_iter.next() {
                    if Some(current_idx) == gt_idx {
                        if val != "." { gt_val = val; }
                    } else if Some(current_idx) == gp_idx {
                        if val != "." {
                            let mut gp_iter = val.split(',');
                            let _gp0 = gp_iter.next();
                            if let Some(v1) = gp_iter.next() { gp1 = v1.parse::<f32>().unwrap_or(0.0); }
                            if let Some(v2) = gp_iter.next() { gp2 = v2.parse::<f32>().unwrap_or(0.0); }
                        }
                    }
                    current_idx += 1;
                }
            }

            let (p0, p1) = if gt_val.starts_with("1|0") {
                (gp2 + gp1, gp2)
            } else if gt_val.starts_with("0|1") {
                (gp2, gp2 + gp1)
            } else {
                (gp2 + (gp1 / 2.0), gp2 + (gp1 / 2.0))
            };

            hap_probs[s][a] = (p0.clamp(1e-5, 1.0 - 1e-5), p1.clamp(1e-5, 1.0 - 1e-5));
        }
    }

    let mut atomic_sample_dists: HashMap<String, Vec<PhasedDist>> = HashMap::new();
    for id in &all_atomic_ids {
        atomic_sample_dists.insert(id.clone(), vec![PhasedDist { p0: 0.0, p1: 0.0 }; num_samples]);
    }

    let mut allele_scores: Vec<(usize, f32)> = Vec::with_capacity(num_alleles);

    for s in 0..num_samples {
        allele_scores.clear();
        for a in 0..num_alleles {
            let score = hap_probs[s][a].0 + hap_probs[s][a].1;
            allele_scores.push((a, score));
        }

        allele_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let m = allele_scores.len().min(max_alleles);
        let top_alleles: Vec<usize> = allele_scores.iter().take(m).map(|x| x.0).collect();

        let mut z0 = 1.0_f32;
        let mut z1 = 1.0_f32;

        let mut w0_map = HashMap::new();
        let mut w1_map = HashMap::new();

        for &a in &top_alleles {
            let (p0, p1) = hap_probs[s][a];
            let w0 = p0 / (1.0 - p0);
            let w1 = p1 / (1.0 - p1);

            w0_map.insert(a, w0);
            w1_map.insert(a, w1);

            z0 += w0;
            z1 += w1;
        }

        for &a in &top_alleles {
            let norm_p0 = w0_map[&a] / z0;
            let norm_p1 = w1_map[&a] / z1;

            for atomic_id in &records[a].atomic_ids {
                if let Some(dists) = atomic_sample_dists.get_mut(atomic_id) {
                    dists[s].p0 += norm_p0;
                    dists[s].p1 += norm_p1;
                }
            }
        }
    }

    let mut sorted_atomic_vars: Vec<_> = all_atomic_ids.into_iter().collect();

    sorted_atomic_vars.sort_by(|a, b| {
        let data_a = id_buffer.get(a);
        let data_b = id_buffer.get(b);

        match (data_a, data_b) {
            (Some(da), Some(db)) => {
                da.0.cmp(&db.0)
                    .then_with(|| da.2.cmp(&db.2))
                    .then_with(|| da.3.cmp(&db.3))
            }
            _ => std::cmp::Ordering::Equal,
        }
    });

    let empty_dists = vec![PhasedDist { p0: 0.0, p1: 0.0 }; num_samples];

    for assigned_id in sorted_atomic_vars {
        let var_data = &id_buffer[&assigned_id];
        let coord = var_data.0;

        let t_rec = records.iter().find(|r| r.atomic_ids.contains(&assigned_id)).unwrap();

        let mut new_info = vec![format!("ID={}", assigned_id)];
        if let Some(v) = &t_rec.raf { new_info.push(format!("RAF={}", v)); }
        if let Some(v) = &t_rec.af { new_info.push(format!("AF={}", v)); }
        if let Some(v) = &t_rec.info_val { new_info.push(format!("INFO={}", v)); }

        let mut vcf_line = vec![
            chrom.clone(),
            coord.to_string(),
            var_data.1.clone(),
            var_data.2.clone(),
            var_data.3.clone(),
            ".".to_string(),
            ".".to_string(),
            new_info.join(";"),
            "GT:DS:GP".to_string(),
        ];

        let dists = atomic_sample_dists.get(&assigned_id).unwrap_or(&empty_dists);

        for s in 0..num_samples {
            let dist = &dists[s];

            let p0 = dist.p0.clamp(0.0, 1.0);
            let p1 = dist.p1.clamp(0.0, 1.0);

            let hap0_gt = if p0 > 0.5 { "1" } else { "0" };
            let hap1_gt = if p1 > 0.5 { "1" } else { "0" };
            let gt = format!("{}|{}", hap0_gt, hap1_gt);

            let ds_str = format_float(p0 + p1);

            let gp0_raw = (1.0 - p0) * (1.0 - p1);
            let gp1_raw = p0 * (1.0 - p1) + (1.0 - p0) * p1;
            let gp2_raw = p0 * p1;

            let mut v0 = (gp0_raw * 1000.0).round() as i32;
            let mut v1 = (gp1_raw * 1000.0).round() as i32;
            let mut v2 = (gp2_raw * 1000.0).round() as i32;

            let diff = 1000 - (v0 + v1 + v2);
            if diff != 0 {
                if v0 >= v1 && v0 >= v2 { v0 += diff; }
                else if v1 >= v0 && v1 >= v2 { v1 += diff; }
                else { v2 += diff; }
            }
            let gp_str = format!("{},{},{}",
                format_float(v0 as f32 / 1000.0),
                format_float(v1 as f32 / 1000.0),
                format_float(v2 as f32 / 1000.0)
            );

            vcf_line.push(format!("{}:{}:{}", gt, ds_str, gp_str));
        }
        writeln!(out_handle, "{}", vcf_line.join("\t")).unwrap();
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: cat <multiallelic VCF> | {} <biallelic ID VCF> <sites VCF> [max_alleles] [window_size]", args[0]);
        std::process::exit(1);
    }

    let vcf_path = &args[1];
    let sites_path = &args[2];
    let max_alleles: usize = if args.len() > 3 {
        args[3].parse().unwrap_or(10)
    } else {
        10
    };

    let window_size: u32 = if args.len() > 4 {
        args[4].parse().unwrap_or(500_000)
    } else {
        500_000
    };

    let mut id_iter = smart_open(vcf_path).lines().peekable();
    let mut sites_iter = smart_open(sites_path).lines();

    let mut next_site = sites_iter.next();
    while let Some(Ok(ref l)) = next_site {
        if l.starts_with('#') {
            next_site = sites_iter.next();
        } else {
            break;
        }
    }

    let stdout = io::stdout();
    let mut out_handle = BufWriter::new(stdout.lock());

    let stdin = io::stdin();
    let mut current_pos: Option<String> = None;
    let mut current_chrom: String = String::new();
    let mut group: Vec<(String, String)> = Vec::new();
    let mut records_processed: usize = 0;

    let mut id_buffer: HashMap<String, (u32, String, String, String)> = HashMap::new();
    let mut active_id_chrom = String::new();

    let start_time = Instant::now();

    eprintln!("Starting Phased Joint-Distribution projection (Max Alleles: {}, Window Size: {})...", max_alleles, window_size);
    for line_result in stdin.lock().lines() {
        let line = line_result.unwrap();
        if line.starts_with('#') {
            if !line.contains("INFO=<ID=AK") && !line.contains("FORMAT=<ID=GL") && !line.contains("FORMAT=<ID=KC") {
                writeln!(out_handle, "{}", line).unwrap();
            }
            continue;
        }

        let site_line = next_site.unwrap_or_else(|| panic!("Error: Sites VCF ran out of records before the main VCF stream")).unwrap();
        next_site = sites_iter.next();

        let site_fields: Vec<&str> = site_line.splitn(9, '\t').collect();
        if site_fields.len() < 5 {
            panic!("Error: Sites VCF line is malformed or missing required columns: {}", site_line);
        }
        let site_chrom = site_fields[0];
        let site_pos = site_fields[1];
        let site_ref = site_fields[3];
        let site_alt = site_fields[4];
        let site_info = if site_fields.len() > 7 { site_fields[7].to_string() } else { "".to_string() };

        let fields: Vec<&str> = line.splitn(6, '\t').collect();
        if fields.len() < 5 { continue; }

        let chrom = fields[0].to_string();
        let pos = fields[1].to_string();
        let ref_seq = fields[3];
        let alt_seq = fields[4];

        if chrom != site_chrom || pos != site_pos || ref_seq != site_ref || alt_seq != site_alt {
            panic!(
                "Error: Lockstep synchronization failed! Variants do not match.\nMain Stream:  {}:{} {} -> {}\nSites Stream: {}:{} {} -> {}",
                chrom, pos, ref_seq, alt_seq, site_chrom, site_pos, site_ref, site_alt
            );
        }

        records_processed += 1;
        if records_processed % 10_000 == 0 {
            let elapsed = start_time.elapsed().as_secs();
            let hours = elapsed / 3600;
            let mins = (elapsed % 3600) / 60;
            let secs = elapsed % 60;
            eprintln!(
                "[{:02}:{:02}:{:02}] Processed {} input records... (Currently at {}:{})",
                hours, mins, secs, records_processed, chrom, pos
            );
        }

        if current_pos.is_none() {
            current_pos = Some(pos.clone());
            current_chrom = chrom.clone();
        }

        if pos != *current_pos.as_ref().unwrap() || chrom != current_chrom {
            process_group(&group, &id_buffer, max_alleles, &mut out_handle);
            group.clear();

            current_pos = Some(pos.clone());
            current_chrom = chrom.clone();
        }

        if group.is_empty() {
            let pos_u32 = pos.parse::<u32>().unwrap_or(0);

            if current_chrom != active_id_chrom {
                id_buffer.clear();
                active_id_chrom = current_chrom.clone();
            }

            id_buffer.retain(|_, v| v.0 + window_size >= pos_u32);

            while let Some(Ok(peek_line)) = id_iter.peek() {
                if peek_line.starts_with('#') {
                    id_iter.next();
                    continue;
                }

                let peek_fields: Vec<&str> = peek_line.splitn(3, '\t').collect();
                let peek_chrom = peek_fields[0];

                if peek_chrom != current_chrom {
                    if active_id_chrom != current_chrom {
                        id_iter.next();
                        continue;
                    } else {
                        break;
                    }
                }

                let peek_pos: u32 = peek_fields[1].parse().unwrap_or(0);

                if peek_pos > pos_u32 + window_size {
                    break;
                }

                if peek_pos + window_size < pos_u32 {
                    id_iter.next();
                    continue;
                }

                let pop_line = id_iter.next().unwrap().unwrap();
                let pop_fields: Vec<&str> = pop_line.split('\t').collect();
                let orig_id = pop_fields[2].to_string();
                let ref_seq = pop_fields[3].to_string();
                let alt_seq = pop_fields[4].to_string();

                for item in pop_fields[7].split(';') {
                    if let Some(id_val) = item.strip_prefix("ID=") {
                        id_buffer.insert(id_val.to_string(), (peek_pos, orig_id, ref_seq, alt_seq));
                        break;
                    }
                }
            }
        }

        group.push((line, site_info));
    }

    if !group.is_empty() {
        process_group(&group, &id_buffer, max_alleles, &mut out_handle);
    }

    out_handle.flush().unwrap();

    let total_elapsed = start_time.elapsed().as_secs();
    let t_hours = total_elapsed / 3600;
    let t_mins = (total_elapsed % 3600) / 60;
    let t_secs = total_elapsed % 60;
    eprintln!(
        "Finished! Processed a total of {} input records in {:02}:{:02}:{:02}.",
        records_processed, t_hours, t_mins, t_secs
    );
}
