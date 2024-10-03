use clap::{Arg, Command};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use bio::io::fastq;

fn main() {
    let matches = Command::new("run_demultiplex_combined")
        .version("1.0")
        .about("Run demultiplexing based on manifest file")
        .arg(
            Arg::new("manifest_file")
                .required(true)
                .index(1)
                .help("Manifest file"),
        )
        .get_matches();

    let manifest_file = matches.get_one::<String>("manifest_file").unwrap();

    let manifest_handle = File::open(manifest_file).expect("Unable to open manifest file");
    let reader = BufReader::new(manifest_handle);

    // Read lines into a Vec and skip the header
    let lines: Vec<String> = reader
        .lines()
        .skip(1) // Skip header line
        .map(|l| l.expect("Error reading manifest file"))
        .collect();

    // Use a thread-safe vector for logging messages
    let log_messages = Arc::new(Mutex::new(Vec::new()));

    // Process each line in parallel
    lines.par_iter().for_each(|line| {
        let line = line.trim();

        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() != 6 {
            eprintln!("Invalid line: {}", line);
            let mut logs = log_messages.lock().unwrap();
            logs.push(format!("Invalid line: {}\n", line));
            return;
        }

        let name = fields[0];
        let file_name = fields[1];
        let _idx1 = fields[2];
        let _seq1 = fields[3];
        let _idx2 = fields[4];
        let seq2 = fields[5];

        // Handle both .fastq and .fastq.gz input files
        let fq_r1_file = if Path::new(&format!("{}_R1_001.fastq.gz", file_name)).is_file() {
            format!("{}_R1_001.fastq.gz", file_name)
        } else if Path::new(&format!("{}_R1_001.fastq", file_name)).is_file() {
            format!("{}_R1_001.fastq", file_name)
        } else {
            eprintln!("R1 file does not exist for {}", file_name);
            let mut logs = log_messages.lock().unwrap();
            logs.push(format!("R1 file does not exist for {}\n", file_name));
            return;
        };

        let fq_r2_file = if Path::new(&format!("{}_R2_001.fastq.gz", file_name)).is_file() {
            format!("{}_R2_001.fastq.gz", file_name)
        } else if Path::new(&format!("{}_R2_001.fastq", file_name)).is_file() {
            format!("{}_R2_001.fastq", file_name)
        } else {
            eprintln!("R2 file does not exist for {}", file_name);
            let mut logs = log_messages.lock().unwrap();
            logs.push(format!("R2 file does not exist for {}\n", file_name));
            return;
        };

        let outbase = format!("{}_{}", name, seq2);

        // Perform demultiplexing directly
        if let Err(e) = demultiplex_fastq_files(&fq_r1_file, &fq_r2_file, seq2, &outbase) {
            eprintln!("Error processing {}: {}", file_name, e);
            let mut logs = log_messages.lock().unwrap();
            logs.push(format!("Error processing {}: {}\n", file_name, e));
        }
    });

    // Write log messages to the logfile
    let logs = log_messages.lock().unwrap();
    if !logs.is_empty() {
        append_to_logfile(&logs.join(""));
    }
}

fn append_to_logfile(message: &str) {
    let mut logfile = OpenOptions::new()
        .create(true)
        .append(true)
        .open("logfile.txt")
        .expect("Unable to open logfile.txt");
    write!(logfile, "{}", message).expect("Unable to write to logfile.txt");
}

fn demultiplex_fastq_files(
    fq_r1_file: &str,
    fq_r2_file: &str,
    adaptseq: &str,
    outbase: &str,
) -> io::Result<()> {
    let adaptseq = adaptseq.as_bytes();
    let index_len = adaptseq.len();

    // Check that files exist
    if !Path::new(fq_r1_file).exists() || !Path::new(fq_r2_file).exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("File(s) do not exist: {}, {}", fq_r1_file, fq_r2_file),
        ));
    }

    // Output files will be compressed (.fastq.gz)
    let outfile1 = format!("{}_L001_R1_001.fastq.gz", outbase);
    let outfile2 = format!("{}_L001_R2_001.fastq.gz", outbase);

    eprintln!("forward read 4N + {}:", String::from_utf8_lossy(adaptseq));
    eprintln!("{} -> {}", fq_r1_file, outfile1);
    eprintln!("{} -> {}", fq_r2_file, outfile2);

    // Open input FASTQ files using bio::io::fastq::Reader
    let in1 = open_fastq_reader(fq_r1_file)?;
    let in2 = open_fastq_reader(fq_r2_file)?;

    // Create output FASTQ writers
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

        // Calculate the start and end indices for slicing
        let start_idx = 4;
        let end_idx = start_idx + index_len;

        // Ensure we have enough sequence length
        if seq1.len() >= end_idx {
            // Perform a case-sensitive comparison
            if &seq1[start_idx..end_idx] == adaptseq {
                count_good += 1;

                // Trim the sequence and quality strings
                let new_seq1 = &seq1[end_idx..];
                let new_qual1 = &qual1[end_idx..];

                // Create new record for R1
                let new_rec1 = fastq::Record::with_attrs(
                    rec1.id(),
                    rec1.desc(),
                    new_seq1,
                    new_qual1,
                );

                // Write to output
                out1.write_record(&new_rec1)?;
                out2.write_record(&rec2)?; // rec2 remains unchanged
            }
        }
    }

    eprintln!("Extracted: {} of {}", count_good, count_total);
    Ok(())
}

fn open_fastq_reader(
    filename: &str,
) -> io::Result<fastq::Reader<Box<dyn BufRead + Send>>> {
    use std::io::Read;

    let file = File::open(filename)?;
    let reader: Box<dyn Read + Send> = if filename.ends_with(".gz") {
        // For gzipped files
        Box::new(MultiGzDecoder::new(file))
    } else {
        // For plain text files
        Box::new(file)
    };

    // Wrap the reader in a BufReader
    let buf_reader = BufReader::new(reader);

    // Box the BufReader as a trait object
    let boxed_buf_reader: Box<dyn BufRead + Send> = Box::new(buf_reader);

    // Use fastq::Reader::from_reader instead of new
    Ok(fastq::Reader::from_bufread(boxed_buf_reader))
}


