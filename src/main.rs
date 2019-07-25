extern crate sufsort_rs;
extern crate bio;
extern crate num;
extern crate clap;

use sufsort_rs::sufsort;
use sufsort_rs::lcp;
use sufsort_rs::rmq;

use bio::io::fasta;
use bio::io::fasta::Record;

use clap::{App, Arg};

const FIRST_DELIMITER: char = '@';
const SECOND_DELIMITER: char = '$';
const FIRST_DELIMITER_U8: u8 = FIRST_DELIMITER as u8;
const SECOND_DELIMITER_U8: u8 = SECOND_DELIMITER as u8;

pub struct RunArgs {
    pub k: usize,
    pub file_name: String,
}

pub struct GST<'s, T> where 
  T: std::marker::Copy + num::Integer + std::fmt::Debug {
    pub first_seq: &'s [u8],
    pub second_seq: &'s [u8],
    pub concat_txt: &'s [u8],
    pub first_eos: usize,
    pub second_eos: usize,
    pub gsa: &'s sufsort::SA<'s, T>,
    pub gisa: &'s Vec<T>,
    pub glcp: &'s Vec<T>,
    pub grmq: &'s rmq::RMQ<'s, T, usize>,
}

pub struct LCPK<T> where 
  T: std::marker::Copy + num::Integer + std::fmt::Debug {
      pub match_length: Vec<T>,
      pub left_match: Vec<T>,
      pub right_match: Vec<T>,
      pub lcpk_length: Vec<T>
}

#[derive(Debug, PartialEq)]
pub enum MatchDirection{
    MatchNone,
    MatchLeft,
    MatchRight
}



pub fn verify_acsk_at<T>(gst: &GST<T>, pos: T, match_pos: T,  max_length: T,
                         k: usize) -> usize
where 
  T: std::marker::Copy + 
     num::Integer + num::FromPrimitive + num::ToPrimitive +
     std::fmt::Debug {
    let mut count: usize = 0;
    let mut upos: usize = T::to_usize(&pos).unwrap();
    let mut umatch_pos: usize = T::to_usize(&match_pos).unwrap();
    let mut max_length: usize = T::to_usize(&max_length).unwrap();
	// // because lcpk = total length - k
    max_length += k;
    while max_length > 0 &&
            upos < gst.concat_txt.len() && umatch_pos < gst.concat_txt.len() {
        if gst.concat_txt[upos] != gst.concat_txt[umatch_pos]{
            count += 1;
        }
        upos += 1; umatch_pos += 1;
        max_length -= 1;
    }
    count
}

pub fn verify_acsk<T>(gst: &GST<T>, lcpk_length: &Vec<T>, match_pos: &Vec<T>,
                      k: usize) -> (usize, usize, usize) where 
  T: std::marker::Copy + 
     num::Integer + num::FromPrimitive + num::ToPrimitive +
     std::fmt::Debug {
	// // prints if # missmatches > k for all i in text 
    let n = gst.concat_txt.len();
    let mut num_less: usize = 0;
    let mut num_more: usize = 0;
    let mut num_equal: usize = 0;
    for i in 2..n {
        let count = verify_acsk_at(gst, gst.gsa.sarray[i], 
                                   gst.gsa.sarray[match_pos[i].to_usize().unwrap()],
                                   lcpk_length[i], k);
        if count > k {
            num_more += 1;
            // println!("{} {} : E# {}", i, gst.gsa.sarray[i].to_usize().unwrap(), count);
        } else if count < k {
            num_less += 1;
            // println!("{} {} : E# {}", i, gst.gsa.sarray[i].to_usize().unwrap(), count);
        } else {
            num_equal += 1;
        }
    }
    println!("No. of matches > {} : {}; No. matches = {}; {} No. matches < {} : {}", 
                k, num_more, k, num_equal, k, num_less);
    (num_more, num_equal, num_less)
}

//
// finding matches for k = 0
//
pub fn find_matches<T>(gst: &GST<T>) -> LCPK<T> where 
  T: std::marker::Copy + num::Integer + num::FromPrimitive + std::fmt::Debug {
    
    let n = gst.concat_txt.len();
    let mut match_length: Vec<T> = vec![T::zero(); n];
    let mut left_match: Vec<T> = vec![T::zero(); n];
    let mut right_match: Vec<T> = vec![T::zero(); n];
    let lcpk_length: Vec<T> = vec![T::zero(); n];
    let eos_first = T::from_usize(gst.first_eos).unwrap();
    let mut min: T = T::zero(); 
    let mut last_match_saidx: T = T::zero();

    // Scan left to right
    for i in 2..(n-1) { // first two entries are delimiters! 
        if (gst.gsa.sarray[i] > eos_first && gst.gsa.sarray[i+1] > eos_first) || 
            (gst.gsa.sarray[i] < eos_first && gst.gsa.sarray[i+1] < eos_first) {
            if gst.glcp[i+1] <= min {
                min = gst.glcp[i+1];
            }
        } else {
            min = gst.glcp[i+1];
            last_match_saidx = T::from_usize(i).unwrap();
        }
        match_length[i+1] = min;
        left_match[i+1] = last_match_saidx;
    }

    min = T::zero(); last_match_saidx = T::zero();
    for i in (3..n).rev()  {
        if (gst.gsa.sarray[i] > eos_first && gst.gsa.sarray[i-1] > eos_first) ||
            (gst.gsa.sarray[i] < eos_first && gst.gsa.sarray[i-1] < eos_first) {
            if gst.glcp[i] <= min {
                min = gst.glcp[i];
            }
        } else {
			min = gst.glcp[i];
			last_match_saidx = T::from_usize(i).unwrap();
        }
        if min > match_length[i-1] {
            match_length[i-1] = min;
            right_match[i-1] = last_match_saidx;
            left_match[i-1] = T::zero();
        } else if min < match_length[i-1] {
            right_match[i-1] = T::zero();
        } else {
            right_match[i-1] = last_match_saidx;
        }
    }

    LCPK::<T>{
        match_length: match_length,
        left_match: left_match,
        right_match: right_match,
        lcpk_length: lcpk_length
    }
}

pub fn forward_match(concat_txt: &[u8], 
                     pos1: usize, pos2: usize,
                     eos_first: usize, eos_second: usize,  
                     k: usize) -> usize {
    let mut lcp: usize = 0;
    let mut rk : i32 = k as i32;
    let mut p1 = pos1;
    let mut p2 = pos2;

    while rk >= 0 && 
            p1 < eos_first && p2 < eos_second &&
            concat_txt[p1] != FIRST_DELIMITER_U8 &&
            concat_txt[p2] != SECOND_DELIMITER_U8  {
        if concat_txt[p1] == concat_txt[p2] {
            lcp += 1; 
        } else {
            rk -= 1;
        }
        p1 += 1; p2 += 1;
    }
    lcp
}

pub fn match_with_second<T>(args: & RunArgs, gst: & GST<T>,
                            lcpk: &LCPK<T>, i: usize) -> (MatchDirection, usize, usize)
    where
    T: std::marker::Copy + 
         num::Integer + num::FromPrimitive + num::ToPrimitive +
        std::fmt::Debug {
    let pos1 = gst.gsa.sarray[i].to_usize().unwrap()  +
                lcpk.match_length[i].to_usize().unwrap();
    let mut pos2: usize = 0;
    let mut maxl: usize = 0;
    let mut flag: MatchDirection = MatchDirection::MatchNone;
    let mut lcp: usize;
    if lcpk.left_match[i] > T::zero() && pos1 < gst.first_eos {
        let mut p = lcpk.left_match[i].to_usize().unwrap();
        while gst.glcp[p+1] >= lcpk.match_length[i] && p > 0 { // TODO: should this be 2?
            if gst.gsa.sarray[p].to_usize().unwrap() > gst.first_eos {
                lcp = forward_match(gst.concat_txt,
                                    pos1,
                                    (gst.gsa.sarray[p] +
                                     lcpk.match_length[i]).to_usize().unwrap(),
                                    gst.first_eos, gst.second_eos, args.k);
                if lcp > maxl {
                    maxl = lcp;
                    pos2 = p;
                    flag = MatchDirection::MatchLeft;
                }
            }
            p -= 1;
        }
    }
    if lcpk.right_match[i] > T::zero() && pos1 < gst.first_eos {
        let mut p = lcpk.right_match[i].to_usize().unwrap();
        while gst.glcp[p] >= lcpk.match_length[i] && p < gst.second_eos {
            if gst.gsa.sarray[p].to_usize().unwrap() > gst.first_eos {
                lcp = forward_match(gst.concat_txt,
                                    pos1,
                                    (gst.gsa.sarray[p] + 
                                        lcpk.match_length[i]).to_usize().unwrap(),
                                    gst.first_eos, gst.second_eos, args.k);
                if lcp > maxl {
                    maxl =lcp;
                    pos2 = p;
                    flag = MatchDirection::MatchRight;
                }
            }
            p += 1;
        }
    }
    (flag, pos2, maxl)
}

pub fn match_with_first<T>(args: & RunArgs, gst: & GST<T>,
                           lcpk: &LCPK<T>, i: usize) -> (MatchDirection, usize, usize)
    where
    T: std::marker::Copy + 
         num::Integer + num::FromPrimitive + num::ToPrimitive +
        std::fmt::Debug {
    let pos1 = gst.gsa.sarray[i].to_usize().unwrap()  +
                lcpk.match_length[i].to_usize().unwrap();

    let mut pos2: usize = gst.second_eos;
    let mut maxl: usize = 0;
    let mut flag: MatchDirection = MatchDirection::MatchNone;
    let mut p: usize; let mut lcp: usize;

    if lcpk.left_match[i] > T::zero() && pos1 < gst.second_eos {
        p = lcpk.left_match[i].to_usize().unwrap();
        while gst.glcp[p+1] >= lcpk.match_length[i] && p > 0 {
            if gst.gsa.sarray[p].to_usize().unwrap() < gst.first_eos {
                lcp = forward_match(gst.concat_txt,
                                    (gst.gsa.sarray[p]+
                                        lcpk.match_length[i]).to_usize().unwrap(), 
                                    pos1,
                                    gst.first_eos, gst.second_eos, args.k);
                if lcp > maxl {
                    maxl = lcp;
                    pos2 = p;
                    flag = MatchDirection::MatchLeft;
                }
            }
            p -= 1;
        }
    }
    if lcpk.right_match[i] > T::zero() && pos1 < gst.second_eos {
        p = lcpk.right_match[i].to_usize().unwrap();
        while gst.glcp[p] >= lcpk.match_length[i] && p < gst.second_eos {
            if gst.gsa.sarray[p].to_usize().unwrap() < gst.first_eos {
                lcp = forward_match(gst.concat_txt,
                                    (gst.gsa.sarray[p] +
                                        lcpk.match_length[i]).to_usize().unwrap(),
                                    pos1,
                                    gst.first_eos, gst.second_eos, args.k);
                if lcp > maxl {
                    maxl = lcp;
                    pos2 = p;
                    flag = MatchDirection::MatchRight;
                }
            }
            p += 1;
        }
    }
    (flag, pos2, maxl)
}

pub fn compute_lcpk_kmacs<T>(args: & RunArgs, gst: & GST<T>, lcpk: &mut LCPK<T>) -> Vec<T>
where 
  T: std::marker::Copy + 
     num::Integer + num::FromPrimitive + num::ToPrimitive +
     std::fmt::Debug {

    assert!(lcpk.lcpk_length.len() == gst.concat_txt.len());
    assert!(lcpk.match_length.len() == gst.concat_txt.len());
    assert!(lcpk.left_match.len() == gst.concat_txt.len());
    assert!(lcpk.right_match.len() == gst.concat_txt.len());

    let mut max_match_kmacs: Vec<T> = vec![T::zero(); gst.gsa.sarray.len()];

    for i in 0..(gst.gsa.sarray.len()) {
        lcpk.lcpk_length[i] = lcpk.match_length[i];

		// forward matching char by char
		let (flag, pos2, maxl) = if gst.gsa.sarray[i].to_usize().unwrap() < gst.first_eos {
            match_with_second(args, gst, lcpk, i)
		} else if gst.gsa.sarray[i].to_usize().unwrap() > gst.first_eos {
            match_with_first(args, gst, lcpk, i)
		} else {
            (MatchDirection::MatchNone, 0, 0)
        };
    
        if flag == MatchDirection::MatchLeft {
            lcpk.left_match[i] = T::from_usize(pos2).unwrap();
            lcpk.right_match[i] = T::zero();
        }
        else if flag == MatchDirection::MatchRight {
            lcpk.left_match[i] = T::zero();
            lcpk.right_match[i] = T::from_usize(pos2).unwrap();
        }
        lcpk.lcpk_length[i] = lcpk.match_length[i] + T::from_usize(maxl).unwrap();

    }
	for i in 0..gst.concat_txt.len() {
		if lcpk.left_match[i] < lcpk.right_match[i] {
			max_match_kmacs[i] = lcpk.right_match[i];
		} else {
			max_match_kmacs[i] = lcpk.left_match[i];
        }
	}
    max_match_kmacs
}

pub fn compute_lcpk_adyar<T>(args: & RunArgs, gst: & GST<T>, lcpk: &mut LCPK<T>,
                            max_match_kmacs: & Vec<T>) -> Vec<T>
where 
  T: std::marker::Copy + 
     num::Integer + num::FromPrimitive + num::ToPrimitive +
     std::fmt::Debug {
    
    assert!(lcpk.lcpk_length.len() == gst.concat_txt.len());
    assert!(lcpk.match_length.len() == gst.concat_txt.len());
    assert!(lcpk.left_match.len() == gst.concat_txt.len());
    assert!(lcpk.right_match.len() == gst.concat_txt.len());
    assert!(max_match_kmacs.len() == gst.concat_txt.len());

    let n = gst.gsa.sarray.len();
    // Phase 1:
	// backward and forward search and storage of missmatches for each i
	let mut max_match_adyar: Vec<T> = max_match_kmacs.clone(); // vec![T::zero(); n];
	for i in 2..n { // skip i = 0 and i = 1 as first two chars are string terminatros
        let mut mismatch_position_bwd: Vec<usize> = vec![n; args.k+1];
        let gsa_pos = gst.gsa.sarray[i].to_usize().unwrap();
        let cur_match_pos = gst.gsa.sarray[max_match_kmacs[i].to_usize().unwrap()].to_usize().unwrap();
		// initialize pos1 , pos2
		let mut pos1 = gsa_pos; 
        let mut pos2 = cur_match_pos;
		let mut kidx = args.k + 1; 

		//backward matching
		while kidx > 0 &&
               ((pos1 < gst.first_eos && pos2 > gst.first_eos) || 
                (pos1 > gst.first_eos && pos2 < gst.first_eos)) {
			if gst.concat_txt[pos1] != gst.concat_txt[pos2] {
				if (pos1 < gst.first_eos && pos2 < gst.first_eos) || 
                    (pos1 > gst.first_eos && pos2 > gst.first_eos) { //checks for '$' crossover
					break;
				}
                mismatch_position_bwd[kidx-1] = pos1 + 1;
				kidx -= 1;
			}
            if pos1 == 0 || pos2 == 0 {
                break;
            }
            pos1 -= 1; pos2 -= 1;
		}

        let mut mismatch_position_fwd: Vec<usize> = vec![n; args.k+1];
		//forward matching
		pos1 = gsa_pos; 
        pos2 = cur_match_pos;
		kidx = args.k + 1;

		while kidx > 0 && pos1 < n && pos2 < n && 
                ((pos1 < gst.first_eos && pos2 > gst.first_eos) || 
                 (pos1 > gst.first_eos && pos2 < gst.first_eos)) {
			if gst.concat_txt[pos1] != gst.concat_txt[pos2] {
				if (pos1 < gst.first_eos && pos2 < gst.first_eos) || 
                    (pos1 > gst.first_eos && pos2 > gst.first_eos) {//checks for '$' crossover
					break;
				}
                mismatch_position_fwd[args.k + 1 - kidx] = pos1 - 1; 
				kidx -= 1;
			}
            pos1 += 1; pos2 += 1;
		}
		// for each missmatch compair sk[i] and [k+j+1]th - [j]th
		// elements of mismatch_position array replace if less

		for j in 0..args.k {
			if mismatch_position_fwd[j] == n || mismatch_position_bwd[j] == n {
				continue;
			}
            assert!(mismatch_position_fwd[j] > mismatch_position_bwd[j]);

            let bw_pos = gst.gisa[mismatch_position_bwd[j]].to_usize().unwrap();
            let lcp_len = lcpk.lcpk_length[bw_pos].to_usize().unwrap();
			//subtract k from substring to comply to with acs definition
            let substring_len = mismatch_position_fwd[j] - mismatch_position_bwd[j] + 1 - args.k;

			if substring_len > lcp_len {
                max_match_adyar[bw_pos] = gst.gisa[mismatch_position_bwd[j]];
                lcpk.lcpk_length[bw_pos] = T::from_usize(substring_len).unwrap();
			}
		}
	}
	// to verify phase 1 of adyar algorithm
	// verify_acsk(gst, &lcpk.lcpk_length, &max_match_adyar, args.k);
    println!( "ACSK PHASE 1 {}", compute_distance(gst, &lcpk));
    //return 0.0;

    // Phase 2:
	// modify sk-array so that sk[i] = x implies sk[i+1] >= x-1
	for i in 1..n {
        let skidx = gst.gisa[gst.gsa.sarray[i].to_usize().unwrap() + 1].to_usize().unwrap();
		if lcpk.lcpk_length[i] > lcpk.lcpk_length[skidx] {
			lcpk.lcpk_length[skidx] = lcpk.lcpk_length[i] - T::one();
			// modify max_match_adyar for phase 2 to account for positions
			// that were not caught by mismatch_postion array.
			max_match_adyar[skidx] = 
                gst.gisa[gst.gsa.sarray[max_match_adyar[i].to_usize().unwrap()].to_usize().unwrap()+1];
		}
	}

    max_match_adyar
}

pub fn compute_distance<T>(gst: &GST<T>, lcpk: &LCPK<T>) -> f64 
    where 
    T: std::marker::Copy + 
        num::Integer + num::FromPrimitive + num::ToPrimitive +
        std::fmt::Debug {
	let mut s1: Vec<usize> = vec![0; gst.first_seq.len()];
    let mut s2: Vec<usize> = vec![0; gst.second_seq.len()];
	for i in 3..gst.gsa.sarray.len() {
        let saidx = gst.gsa.sarray[i].to_usize().unwrap();
		if saidx > gst.first_eos {
			s2[saidx-gst.first_eos] = lcpk.lcpk_length[i].to_usize().unwrap();
		}
        if saidx < gst.first_seq.len() {
			s1[saidx] = lcpk.lcpk_length[i].to_usize().unwrap();
	    }
    }

	// calculate avgs1 and svgs2 and d_acs
	let mut avg_s1: f64=0.0; let mut avg_s2: f64=0.0;
	for i in 0..(s1.len()) {avg_s1 += s1[i] as f64;}
	avg_s1 = avg_s1/(s1.len()) as f64;
	for i in 0..(s2.len()) {avg_s2 += s2[i] as f64;}
	avg_s2 = avg_s2/s2.len() as f64;

	let d_acs: f64 =
     ((( (s1.len()) as f64).log10()/(2.0*avg_s2)) + ((s2.len() as f64).log10()/(2.0*avg_s1))) - 
        ((((s1.len()) as f64).log10()/(((s1.len())  as f64))) + ((s2.len() as f64).log10()/(s2.len() as f64)));
	//std::cout << "ACS for k= " << args.k << " : " << d_acs << std::endl;
	// args.d_vector.push_back(d_acs);
    d_acs
}


pub fn compute_acsk<T>(args: & RunArgs, gst: & GST<T>, lcpk: &mut LCPK<T>) -> f64
    where 
    T: std::marker::Copy + 
        num::Integer + num::FromPrimitive + num::ToPrimitive +
        std::fmt::Debug {

	// lcpk based on k-macs's heuristic
	// to calculate and print acsk based in k-macs
	let max_match_kmacs: Vec<T> = compute_lcpk_kmacs(args, gst, lcpk);


	// to verify k-macs algorithm
	//verify_acsk(gst, &lcpk.lcpk_length, &max_match, args.k);
    let _acs_kmacs = compute_distance(gst, &lcpk);

    // run adyar algorithm using kmacs inputs
    let _max_match_adyar = compute_lcpk_adyar(args, gst, lcpk, &max_match_kmacs);

	// to verify adyar algorithm
	// verify_acsk(gst, &lcpk.lcpk_length, &_max_match_adyar, args.k);
    compute_distance(gst, &lcpk)
}


fn run_adyar(rargs: RunArgs) {
    let reader = fasta::Reader::from_file(rargs.file_name.as_str()).unwrap();
    let rcds: Vec<Record> = reader.records().map(|x| x.unwrap()).collect();
    println!("No. of Sequences {}", rcds.len());
    let nseq = rcds.len();

    for x in 0..nseq{
        for y in (x+1)..nseq {
            let seqx = rcds[x].seq();
            let seqy = rcds[y].seq();
            let mut tvec: Vec<u8> = seqx.to_vec();
            let eos_first = tvec.len();
            tvec.push(FIRST_DELIMITER_U8);
            tvec.append(&mut seqy.to_vec());
            let eos_second = tvec.len();
            tvec.push(SECOND_DELIMITER_U8);
            // let tlen = tvec.len();

            let gsa = sufsort::SA::<i32>::new(&tvec);
            let gisa = sufsort::construct_isa(&gsa.sarray);
            let glcp = lcp::construct_lcp_from_sa(gsa.txt, &gsa.sarray, &gisa);
            let grmq = rmq::RMQ::<i32, usize>::new(&glcp);

            let gst: GST<i32> = GST::<i32> {    
                first_seq: &seqx, second_seq: &seqy, concat_txt: &tvec,
                first_eos: eos_first, second_eos: eos_second,
                gsa: &gsa, gisa: &gisa, glcp: &glcp, grmq: &grmq
            };

            println!("First Len {}, Second Len {}, First EOS {}, Second EOS {}", 
                        gst.first_seq.len(), gst.second_seq.len(),
                        gst.first_eos, gst.second_eos);

            let mut lcpk = find_matches(&gst);
            // println!("Matches {:?}", lcpk.match_length);
            let acsk = compute_acsk(&rargs, &gst, &mut lcpk);

            println!("{}", acsk);
        }
    // }
    }

}

fn main() {
    fn is_integer(v: String) -> Result<(), String> {
        if v.parse::<usize>().is_ok() { return Ok(()); }
        Err(String::from("The value is not valid integer"))
    }
    let matches = App::new("Adyar : Heuristic")
                            .version("0.1.0")
                            .author("Sriram P C")
                            .about("Does awesome things")
                            .arg(Arg::with_name("INPUT")
                                .help("Sets the input file to use")
                                .required(true)
                                .index(1))
                            .arg(Arg::with_name("mismatches")
                               .short("k")
                               .long("mismatches")
                               .value_name("NO_MISMATCHES")
                               .validator(is_integer)
                               .default_value("1")
                               .help("No. of Mismatches")
                               .takes_value(true))
                          .get_matches();

    
    let rargs: RunArgs = RunArgs {
        file_name: String::from(matches.value_of("INPUT").unwrap()),
        k: matches.value_of("mismatches").unwrap().parse::<usize>().unwrap(),
    };
    println!("Using input args: INPUT = {}, K = {} ", rargs.file_name, rargs.k);
    run_adyar(rargs);
}
