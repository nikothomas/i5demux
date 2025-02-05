use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

use bio::io::fastq;
use clap::Parser;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

/// A simple CLI for running demultiplexing based on a manifest file
#[derive(Parser, Debug)]
#[command(
    name = "run_demultiplex_combined",
    version = "1.0",
    about = "Run demultiplexing based on manifest file"
)]
struct Cli {
    /// Path to the manifest file (tab-delimited with 6 columns; the first line is a header)
    manifest_file: String,
}

fn main() -> io::Result<()> {
    // Parse CLI arguments
    let args = Cli::parse();

    // Read manifest lines
    let manifest_file = &args.manifest_file;
    let manifest_handle = File::open(manifest_file)
        .unwrap_or_else(|_| panic!("Unable to open manifest file: {}", manifest_file));
    let reader = BufReader::new(manifest_handle);

    // Collect lines (skip header)
    let lines: Vec<String> = reader
        .lines()
        .skip(1) // Skip header line
        .map(|l| l.expect("Error reading manifest file"))
        .collect();

    // Setup a thread-safe logger (for everything that was previously eprinted)
    let log_messages = Arc::new(Mutex::new(Vec::new()));

    // Create a progress bar (just for progress display; no printing)
    let pb = Arc::new(
        ProgressBar::new(lines.len() as u64)
            .with_message("Processing manifest lines..."),
    );

    // Customize the progress bar style
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>3}/{len:3} {msg}")
            .unwrap()
            .progress_chars("=> "),
    );

    // Process each line in parallel
    lines.par_iter().for_each(|line| {
        let pb = Arc::clone(&pb);
        let log_messages = Arc::clone(&log_messages);

        let line = line.trim();
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() != 6 {
            let mut logs = log_messages.lock().unwrap();
            logs.push(format!("Invalid line: {}\n", line));
            pb.inc(1);
            return;
        }

        let name = fields[0];
        let file_name = fields[1];
        let _idx1 = fields[2];
        let _seq1 = fields[3];
        let _idx2 = fields[4];
        let seq2 = fields[5];

        // Figure out the file names
        let fq_r1_file = if Path::new(&format!("{}_R1_001.fastq.gz", file_name)).is_file() {
            format!("{}_R1_001.fastq.gz", file_name)
        } else if Path::new(&format!("{}_R1_001.fastq", file_name)).is_file() {
            format!("{}_R1_001.fastq", file_name)
        } else {
            let mut logs = log_messages.lock().unwrap();
            logs.push(format!("R1 file does not exist for {}\n", file_name));
            pb.inc(1);
            return;
        };

        let fq_r2_file = if Path::new(&format!("{}_R2_001.fastq.gz", file_name)).is_file() {
            format!("{}_R2_001.fastq.gz", file_name)
        } else if Path::new(&format!("{}_R2_001.fastq", file_name)).is_file() {
            format!("{}_R2_001.fastq", file_name)
        } else {
            let mut logs = log_messages.lock().unwrap();
            logs.push(format!("R2 file does not exist for {}\n", file_name));
            pb.inc(1);
            return;
        };

        let outbase = format!("{}_{}", name, seq2);

        // Perform demultiplexing
        if let Err(e) = demultiplex_fastq_files(
            &fq_r1_file,
            &fq_r2_file,
            seq2,
            &outbase,
            &log_messages
        ) {
            let mut logs = log_messages.lock().unwrap();
            logs.push(format!("Error processing {}: {}\n", file_name, e));
        }

        // Once done with this line, increment the progress bar
        pb.inc(1);
    });

    // Finalize the progress bar
    pb.finish_with_message("Done processing lines");

    // Write log messages to logfile
    let logs = log_messages.lock().unwrap();
    if !logs.is_empty() {
        append_to_logfile(&logs.join(""));
    }

    Ok(())
}

/// Appends `message` to `logfile.txt`.
fn append_to_logfile(message: &str) {
    let mut logfile = OpenOptions::new()
        .create(true)
        .append(true)
        .open("logfile.txt")
        .expect("Unable to open logfile.txt");
    write!(logfile, "{}", message).expect("Unable to write to logfile.txt");
}

/// Demultiplex FASTQ files without printing anything to stdout/stderr; logs go to `log_messages`.
fn demultiplex_fastq_files(
    fq_r1_file: &str,
    fq_r2_file: &str,
    adaptseq: &str,
    outbase: &str,
    log_messages: &Arc<Mutex<Vec<String>>>,
) -> io::Result<()> {
    let adaptseq = adaptseq.as_bytes();
    let index_len = adaptseq.len();

    // Check file existence
    if !Path::new(fq_r1_file).exists() || !Path::new(fq_r2_file).exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("File(s) do not exist: {}, {}", fq_r1_file, fq_r2_file),
        ));
    }

    // Prepare output file names
    let outfile1 = format!("{}_L001_R1_001.fastq.gz", outbase);
    let outfile2 = format!("{}_L001_R2_001.fastq.gz", outbase);

    // Log the demultiplexing action
    {
        let mut logs = log_messages.lock().unwrap();
        logs.push(format!(
            "Demultiplexing:\n  {} -> {}\n  {} -> {}\n",
            fq_r1_file, outfile1, fq_r2_file, outfile2
        ));
    }

    let in1 = open_fastq_reader(fq_r1_file)?;
    let in2 = open_fastq_reader(fq_r2_file)?;

    let out1_file = File::create(&outfile1)?;
    let out2_file = File::create(&outfile2)?;

    let encoder1 = GzEncoder::new(out1_file, Compression::default());
    let encoder2 = GzEncoder::new(out2_file, Compression::default());

    let mut out1 = fastq::Writer::new(encoder1);
    let mut out2 = fastq::Writer::new(encoder2);

    let mut records1 = in1.records();
    let mut records2 = in2.records();

    let mut count_good = 0;
    let mut count_total = 0;

    loop {
        let rec1 = match records1.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(io::Error::new(io::ErrorKind::Other, e)),
            None => break,
        };
        let rec2 = match records2.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(io::Error::new(io::ErrorKind::Other, e)),
            None => break,
        };

        count_total += 1;

        let seq1 = rec1.seq();
        let qual1 = rec1.qual();

        let start_idx = 4;
        let end_idx = start_idx + index_len;

        // If the read in R1 matches the adapter after the first 4 bases, trim
        if seq1.len() >= end_idx && &seq1[start_idx..end_idx] == adaptseq {
            count_good += 1;

            let new_seq1 = &seq1[end_idx..];
            let new_qual1 = &qual1[end_idx..];
            let new_rec1 = fastq::Record::with_attrs(rec1.id(), rec1.desc(), new_seq1, new_qual1);

            out1.write_record(&new_rec1)?;
            out2.write_record(&rec2)?;
        }
    }

    {
        let mut logs = log_messages.lock().unwrap();
        logs.push(format!("Extracted: {} of {}\n", count_good, count_total));
    }

    Ok(())
}

/// Opens a FASTQ file (gzipped or not) as a `bio::io::fastq::Reader`.
fn open_fastq_reader(filename: &str) -> io::Result<fastq::Reader<Box<dyn BufRead + Send>>> {
    use std::io::Read;

    let file = File::open(filename)?;
    let reader: Box<dyn Read + Send> = if filename.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let buf_reader = BufReader::new(reader);
    Ok(fastq::Reader::from_bufread(Box::new(buf_reader)))
}
