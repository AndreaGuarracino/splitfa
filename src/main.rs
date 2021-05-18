#![allow(clippy::too_many_arguments)]

use std::io::{self};
use seq_io::fasta::{Reader, Record};
//use substring::Substring;
use std::str;

extern crate clap;

use clap::{App, Arg};

use rand_distr::{Beta, Distribution, BetaError};
use rand::thread_rng;

fn split_fasta(input: &str, seg_length_min: usize, seg_length_max: usize, step: f32) -> Result<(), BetaError>{
    let seg_length_range = seg_length_max - seg_length_min;
    let mut rng = thread_rng();
    let beta = Beta::new(1.5, 15.0).unwrap();
    let mut reader = Reader::from_path(input).unwrap();
    while let Some(result) = reader.next() {
        let record = result.unwrap();
        let seq = record.full_seq();
        let name = record.id().unwrap();
        let mut start: usize = 0;
        let total_length: usize = seq.len();
        let mut seg_length = seg_length_min + (beta.sample(&mut rng) * seg_length_range as f64) as usize;
        if total_length < seg_length {
            println!(">{}:{}-{}", name, 0, total_length);
            println!("{}", str::from_utf8(&seq[0..total_length]).unwrap());
        } else {
            while start + seg_length <= total_length {
                println!(">{}:{}-{}", name, start, start + seg_length);
                println!("{}", str::from_utf8(&seq[start..(start + seg_length)]).unwrap());
                start += (step * seg_length as f32) as usize;

                seg_length = seg_length_min + (beta.sample(&mut rng) * seg_length_range as f64) as usize;
            }
            if start < total_length {
                start = total_length - seg_length;
                println!(">{}:{}-{}", name, start, start + seg_length);
                println!("{}", str::from_utf8(&seq[start..(start + seg_length)]).unwrap());
            }
        }
    }

    Ok(())
}

fn main() -> io::Result<()> {
    let matches = App::new("splitfa")
        .version("0.1.0")
        .author("Erik Garrison <erik.garrison@gmail.com>")
        .about("Split a FASTA file into subsequences of a given length and overlap length")
        .arg(
            Arg::with_name("INPUT")
                .required(true)
                .takes_value(true)
                .index(1)
                .help("input FASTA file"),
        )
        .arg(
            Arg::with_name("seg-length")
                .short("l")
                .long("seg-length")
                .takes_value(true)
                .help("Random length of the splits, following a Poisson distribution: mean:min-max"),
        )
        .arg(
            Arg::with_name("step")
                .short("s")
                .long("step")
                .takes_value(true)
                .help("Step size between splits"),
        )
        .get_matches();

    let filename = matches.value_of("INPUT").unwrap();

    let seg_length_min = matches.value_of("seg-length")
        .unwrap()
        .split('-')
        .next()
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let seg_length_max = matches.value_of("seg-length")
        .unwrap()
        .split('-')
        .nth(1)
        .unwrap()
        .parse::<usize>()
        .unwrap();

    let step = matches.value_of("step").unwrap().parse::<f32>().unwrap();

    split_fasta(filename, seg_length_min, seg_length_max, step);

    Ok(())
}
