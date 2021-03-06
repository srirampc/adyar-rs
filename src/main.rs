extern crate sufsort_rs;
extern crate bio;
extern crate num;
extern crate clap;

use std::io::prelude::*;
use std::fs::File;

use sufsort_rs::sufsort;
use sufsort_rs::lcp;
use sufsort_rs::rmq;

use bio::io::fasta;
use bio::io::fasta::Record;

use clap::{App, Arg};

const FIRST_DELIMITER: char = '#';
const SECOND_DELIMITER: char = '$';
const FIRST_DELIMITER_U8: u8 = FIRST_DELIMITER as u8;
const SECOND_DELIMITER_U8: u8 = SECOND_DELIMITER as u8;

pub struct RunArgs {
    pub k: usize,
    pub file_name: String,
    pub out_file: String,
    pub use_rmq: bool,
    pub show_progress: bool,
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
    pub grmq: &'s Option<rmq::RMQ<'s, T, usize>>,
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



pub fn small_hist<T>(invec: & Vec<T>)
    where T: num::Integer + num::FromPrimitive + num::ToPrimitive{
    let hmax = 8;
    let mut shist : Vec<usize> = vec![0; hmax];
    for i in 0..invec.len() {
        if invec[i].to_usize().unwrap() < hmax {
            shist[invec[i].to_usize().unwrap()] += 1;
        } 
    }

    println!("Small Histogram");
    for i in 0..shist.len() {
        println!("{} : {}", i, shist[i]);
    }
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
pub fn find_exact_matches<T>(gst: &GST<T>) -> LCPK<T> where 
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
    for i in 0..(n-1) { // first two entries are delimiters! 
        if (gst.gsa.sarray[i] > eos_first && gst.gsa.sarray[i+1] > eos_first) || 
            (gst.gsa.sarray[i] < eos_first && gst.gsa.sarray[i+1] < eos_first) {
            if gst.glcp[i+1] < min {
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
    for i in (1..n).rev()  {
        if (gst.gsa.sarray[i] > eos_first && gst.gsa.sarray[i-1] > eos_first) ||
            (gst.gsa.sarray[i] < eos_first && gst.gsa.sarray[i-1] < eos_first) {
            if gst.glcp[i] < min {
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
                     k: usize) -> (usize, usize) {
    let mut lcp: usize = 0;
    let mut lcpkm1: usize = 0;
    let mut rk : i32 = k as i32;
    let mut p1 = pos1;
    let mut p2 = pos2;

    while rk >= 0 && 
            p1 < eos_first && p2 < eos_second  {
        if concat_txt[p1] == concat_txt[p2] {
            lcp += 1; 
        } else {
            rk -= 1;
            if rk == 0 {
                lcpkm1 = lcp;
            }
        }
        p1 += 1; p2 += 1;
    }
    (lcp, lcpkm1)
}

pub fn match_with_second<T>(args: & RunArgs, gst: & GST<T>,
                            lcpk: &LCPK<T>, i: usize) -> 
                (MatchDirection, usize, usize,
                 MatchDirection, usize, usize, usize)
    where
    T: std::marker::Copy + 
         num::Integer + num::FromPrimitive + num::ToPrimitive +
        std::fmt::Debug {
    let pos1 = gst.gsa.sarray[i].to_usize().unwrap()  +
                lcpk.match_length[i].to_usize().unwrap();
    let mut pos2: usize = 0;
    let mut maxl: usize = 0;
    let mut flag: MatchDirection = MatchDirection::MatchNone;
    let mut pos2km: usize = 0;
    let mut maxlkm: usize = 0;
    let mut flagkm: MatchDirection = MatchDirection::MatchNone;
    let mut maxp: usize = 0;
    if lcpk.left_match[i] > T::zero() && pos1 < gst.first_eos {
        let mut p = lcpk.left_match[i].to_usize().unwrap();
        while gst.glcp[p+1] >= lcpk.match_length[i] && p > 0 { // TODO: should this be 2?
            if gst.gsa.sarray[p].to_usize().unwrap() > gst.first_eos {
                let (lcp, lcpkm) = forward_match(gst.concat_txt,
                                            pos1,
                                            (gst.gsa.sarray[p] +
                                            lcpk.match_length[i]).to_usize().unwrap(),
                                            gst.first_eos, gst.second_eos, args.k);
                if lcp > maxl {
                    maxl = lcp;
                    pos2 = p;
                    flag = MatchDirection::MatchLeft;
                }
                if lcpkm > maxlkm {
                    maxlkm = lcpkm;
                    pos2km = p;
                    flagkm = MatchDirection::MatchLeft;
                }
            }
            p -= 1;
        }
        maxp = lcpk.left_match[i].to_usize().unwrap() - p;
    }
    if lcpk.right_match[i] > T::zero() && pos1 < gst.first_eos {
        let mut p = lcpk.right_match[i].to_usize().unwrap();
        while gst.glcp[p] >= lcpk.match_length[i] && p < gst.second_eos {
            if gst.gsa.sarray[p].to_usize().unwrap() > gst.first_eos {
                let (lcp, lcpkm) = forward_match(gst.concat_txt,
                                    pos1,
                                    (gst.gsa.sarray[p] + 
                                        lcpk.match_length[i]).to_usize().unwrap(),
                                    gst.first_eos, gst.second_eos, args.k);
                if lcp > maxl {
                    maxl =lcp;
                    pos2 = p;
                    flag = MatchDirection::MatchRight;
                }
                if lcpkm > maxlkm {
                    maxlkm = lcpkm;
                    pos2km = p;
                    flagkm = MatchDirection::MatchRight;
                }
            }
            p += 1;
        }
        maxp = std::cmp::max(maxp, p - lcpk.right_match[i].to_usize().unwrap());
    }
    (flag, pos2, maxl, flagkm, pos2km, maxlkm, maxp)
}

pub fn match_with_first<T>(args: & RunArgs, gst: & GST<T>,
                           lcpk: &LCPK<T>, i: usize) ->
                (MatchDirection, usize, usize,
                 MatchDirection, usize, usize, usize)
    where
    T: std::marker::Copy + 
         num::Integer + num::FromPrimitive + num::ToPrimitive +
        std::fmt::Debug {
    let pos1 = gst.gsa.sarray[i].to_usize().unwrap()  +
                lcpk.match_length[i].to_usize().unwrap();

    let mut pos2: usize = gst.second_eos;
    let mut maxl: usize = 0;
    let mut flag: MatchDirection = MatchDirection::MatchNone;
    let mut pos2km: usize = 0;
    let mut maxlkm: usize = 0;
    let mut flagkm: MatchDirection = MatchDirection::MatchNone;
    let mut maxp: usize = 0; 

    if lcpk.left_match[i] > T::zero() && pos1 < gst.second_eos {
        let mut p = lcpk.left_match[i].to_usize().unwrap();
        while gst.glcp[p+1] >= lcpk.match_length[i] && p > 0 {
            if gst.gsa.sarray[p].to_usize().unwrap() < gst.first_eos {
                let (lcp, lcpkm) = forward_match(gst.concat_txt,
                                    (gst.gsa.sarray[p]+
                                        lcpk.match_length[i]).to_usize().unwrap(), 
                                    pos1,
                                    gst.first_eos, gst.second_eos, args.k);
                if lcp > maxl {
                    maxl = lcp;
                    pos2 = p;
                    flag = MatchDirection::MatchLeft;
                }
                if lcpkm > maxlkm {
                    maxlkm = lcpkm;
                    pos2km = p;
                    flagkm = MatchDirection::MatchLeft;
                }
            }
            p -= 1;
        }
        maxp = lcpk.left_match[i].to_usize().unwrap() - p;
    }
    if lcpk.right_match[i] > T::zero() && pos1 < gst.second_eos {
        let mut p = lcpk.right_match[i].to_usize().unwrap();
        while gst.glcp[p] >= lcpk.match_length[i] && p < gst.second_eos {
            if gst.gsa.sarray[p].to_usize().unwrap() < gst.first_eos {
                let (lcp, lcpkm) = forward_match(gst.concat_txt,
                                    (gst.gsa.sarray[p] +
                                        lcpk.match_length[i]).to_usize().unwrap(),
                                    pos1,
                                    gst.first_eos, gst.second_eos, args.k);
                if lcp > maxl {
                    maxl = lcp;
                    pos2 = p;
                    flag = MatchDirection::MatchRight;
                }
                if lcpkm > maxlkm {
                    maxlkm = lcpkm;
                    pos2km = p;
                    flagkm = MatchDirection::MatchRight;
                }
            }
            p += 1;
        }
        maxp = std::cmp::max(maxp, p - lcpk.right_match[i].to_usize().unwrap());
    }
    (flag, pos2, maxl, flagkm, pos2km, maxlkm, maxp)
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
    let mut loop_maxp = 0;

    for i in 2..(gst.gsa.sarray.len()) {
        lcpk.lcpk_length[i] = lcpk.match_length[i];
        if lcpk.match_length[i] == T::zero() {
            continue;
        }
        let gpos = gst.gsa.sarray[i].to_usize().unwrap();

		// forward matching char by char
        assert!(gpos != gst.first_eos && gpos != gst.second_eos);
		let (flag, pos2, maxl,
             flagkm, pos2km, maxlkm, maxp) =
            if gpos < gst.first_eos {
                match_with_second(args, gst, lcpk, i)
            } else {
                match_with_first(args, gst, lcpk, i)
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
        loop_maxp = std::cmp::max(loop_maxp, maxp);

        if gpos > 1 && gpos != gst.first_eos + 1 && gpos != gst.second_eos + 1  {
            let pgpos = gpos - 1;
            let pgidx = gst.gisa[pgpos].to_usize().unwrap();
            if pgidx > 2 && pgidx <  gst.concat_txt.len() &&
                 lcpk.match_length[pgidx] == T::zero() {
                if flagkm == MatchDirection::MatchLeft {
                    lcpk.left_match[pgidx] = T::from_usize(pos2km).unwrap();
                    lcpk.right_match[pgidx] = T::zero();
                }
                else if flag == MatchDirection::MatchRight {
                    lcpk.left_match[pgidx] = T::zero();
                    lcpk.right_match[pgidx] = T::from_usize(pos2km).unwrap();
                }
                lcpk.lcpk_length[pgidx] = lcpk.match_length[pgidx] + T::from_usize(maxlkm).unwrap();
            }   
        }
    }
    // println!("Loop maxp {}", loop_maxp);

    for i in (0..gst.gsa.sarray.len()).rev(){
        if i == gst.first_eos || i == gst.second_eos {
            continue;
        }
        let x = gst.gisa[i].to_usize().unwrap();
        if i+1 == gst.first_eos || i+1 == gst.second_eos {
            // TODO
            continue;
        }
        assert!(x >= 2);
        if lcpk.match_length[x] == T::zero() && lcpk.left_match[i] == T::zero() &&
             lcpk.left_match[x] == lcpk.right_match[x] {
            let y = gst.gisa[i + 1].to_usize().unwrap();
            if lcpk.left_match[y] != T::zero() {
                let (ml, _a) = forward_match(gst.concat_txt, i+1,
                                       lcpk.left_match[y].to_usize().unwrap() +1, 
                                        gst.first_eos, 
                                        gst.second_eos, args.k - 1);
                lcpk.match_length[x] = T::from_usize(ml).unwrap();
                lcpk.left_match[x] = T::from_usize(y).unwrap();
            } else if lcpk.right_match[y] != T::zero(){
                let (ml, _a) = forward_match(gst.concat_txt, i+1, 
                                        lcpk.right_match[y].to_usize().unwrap()+1, 
                                        gst.first_eos, 
                                        gst.second_eos, args.k - 1);
                lcpk.match_length[x] = T::from_usize(ml).unwrap();
                lcpk.right_match[x] = T::from_usize(y).unwrap();
            }
        }
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
            //assert!(mismatch_position_fwd[j] > mismatch_position_bwd[j]);
            if mismatch_position_fwd[j] < mismatch_position_bwd[j] + args.k { continue };

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
    if args.show_progress == true {
        print!("1");
    }

    #[cfg(debug_assertions)]
    println!( "ACSK PHASE 1 {}", compute_distance(gst, &lcpk));
    //return 0.0;

    // Phase 2:
	// modify sk-array so that sk[i] = x implies sk[i+1] >= x-1
	for i in 2..n {
        let skidx = gst.gisa[gst.gsa.sarray[i].to_usize().unwrap() + 1].to_usize().unwrap();
		if lcpk.lcpk_length[i] > lcpk.lcpk_length[skidx] {
			lcpk.lcpk_length[skidx] = lcpk.lcpk_length[i] - T::one();
			// modify max_match_adyar for phase 2 to account for positions
			// that were not caught by mismatch_postion array.
			max_match_adyar[skidx] = 
                gst.gisa[gst.gsa.sarray[max_match_adyar[i].to_usize().unwrap()].to_usize().unwrap()+1];
		}
	}

    if args.show_progress == true {
        print!("2");
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
	for i in 2..gst.gsa.sarray.len() {
        let saidx = gst.gsa.sarray[i].to_usize().unwrap();
		if saidx > gst.first_eos {
			s2[saidx-gst.first_eos-1] = lcpk.lcpk_length[i].to_usize().unwrap();
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
     (
         ((s1.len() as f64).log10()/(2.0 * avg_s2)) + ((s2.len() as f64).log10()/(2.0*avg_s1))
        ) - 
     (
         ((s1.len() as f64).log10()/(s1.len()  as f64)) + ((s2.len() as f64).log10()/(s2.len() as f64))
        );
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
    
    if args.show_progress == true {
        print!("K"); // indicate kmacs part is done
    }


	// to verify k-macs algorithm
	//verify_acsk(gst, &lcpk.lcpk_length, &max_match, args.k);
    #[cfg(debug_assertions)]
    let _acs_kmacs = compute_distance(gst, &lcpk);

    #[cfg(debug_assertions)]
    println!( "ACSK KMACS {}", compute_distance(gst, &lcpk));

    // run adyar algorithm using kmacs inputs
    let _max_match_adyar = compute_lcpk_adyar(args, gst, lcpk, &max_match_kmacs);

	// to verify adyar algorithm
	// verify_acsk(gst, &lcpk.lcpk_length, &_max_match_adyar, args.k);
    let rdist = compute_distance(gst, &lcpk);

    if args.show_progress == true {
        print!("C");
    }
    #[cfg(debug_assertions)]
    println!( "ACSK ADYAR {}", rdist);

    rdist
}


pub fn fill_dist_matrix(nseq: usize, pairwise_dcs:& Vec<f64>) -> Vec<Vec<f64>> {
	let mut d_matrix: Vec<Vec<f64>> = vec![vec![0.0; nseq]; nseq];
	//populate d_matrix using d_vector
	let mut idx: usize = 0;

	for i in 0..nseq {
		for j in (i+1)..nseq {
			let value = pairwise_dcs[idx];
			//std::cout << value << std::endl;
			d_matrix[i][j] = value;
			d_matrix[j][i] = value;
			idx += 1;
		}
	}
    d_matrix
}

pub fn write_dist_matrix(args : &RunArgs,
                         d_matrix: & Vec<Vec<f64>>,
                         title: &Vec<String>) -> std::io::Result<()> {
    let mut outf = File::create(&args.out_file)?;
    //let mut writer = BufWriter::new(outf);
    let nseq: usize = d_matrix.len();
    let mut mat_string = String::new();
    mat_string.push_str(&nseq.to_string());
    mat_string.push_str("\n");
	//file_object << nReads << std::endl;

	for i in 0..nseq {
		let name = &title[i];
		// if(name.len() > 10)
		// 	name = name.(0,10);
        let nlen = if name.len() > 10 {
            let mut txstr = name.to_string();
            txstr.truncate(10);
            mat_string.push_str(&txstr);
            txstr.len()
        } else {
            mat_string.push_str(name);
            name.len()
        };
        for _j in 0..(14-nlen) {
            mat_string.push(' ');
        }
		// file_object << std::setw(14) << std::left << name;
        for j in 0..nseq {
			let kdxy = d_matrix[i][j];
			// file_object << kdxy << ((j==nReads-1)?"":" ");
            mat_string.push_str(&kdxy.to_string());
            if j < nseq - 1 {
                mat_string.push(' ');
            }
		}
		// file_object << std::endl;
        mat_string.push('\n');
	}

    outf.write_all(mat_string.as_bytes())?;
    Ok(())
}

pub fn write_output(args: &RunArgs, titles: &Vec<String>,
                    pairwise_dcs: &Vec<f64>) -> std::io::Result<()> {
    let d_matrix = fill_dist_matrix(titles.len(), pairwise_dcs);
    write_dist_matrix(args, &d_matrix, titles)?;
    Ok(())
}

pub fn adyar_pairwise_acsk(rargs: &RunArgs, seqx: & [u8], seqy: & [u8]) -> f64 {
    let mut tvec: Vec<u8> = seqx.to_vec();
    let eos_first = tvec.len();
    tvec.push(FIRST_DELIMITER_U8);
    let mut yvec: Vec<u8> = seqy.to_vec();
    tvec.append(&mut yvec);
    let eos_second = tvec.len();
    tvec.push(SECOND_DELIMITER_U8);
    let gsa = sufsort::SA::<i32>::new(&tvec);
    let gisa = sufsort::construct_isa(&gsa.sarray);
    let glcp = lcp::construct_lcp_from_sa(gsa.txt, &gsa.sarray, &gisa);
    let grmq = if rargs.use_rmq {
            Some(rmq::RMQ::<i32, usize>::new(&glcp))
        } else {
            None
        };

    #[cfg(debug_assertions)]
    small_hist(&glcp);

    #[cfg(debug_assertions)]
    println!("Construction Complete {}", tvec.len());

    let gst = GST::<i32> {    
        first_seq: &seqx, second_seq: &seqy, concat_txt: &tvec,
        first_eos: eos_first, second_eos: eos_second,
        gsa: &gsa, gisa: &gisa, glcp: &glcp, grmq: &grmq
    };

    // println!("First Len {}, Second Len {}, First EOS {}, Second EOS {}", 
    //             gst.first_seq.len(), gst.second_seq.len(),
    //             gst.first_eos, gst.second_eos);

    let mut lcpk = find_exact_matches(&gst);

    #[cfg(debug_assertions)]
    small_hist(&lcpk.match_length);

    let acsk = compute_acsk(&rargs, &gst, &mut lcpk);
    acsk
}


fn run_adyar(rargs: RunArgs) -> std::io::Result<()> {

    let reader = fasta::Reader::from_file(rargs.file_name.as_str()).unwrap();
    let rcds: Vec<Record> = reader.records().map(|x| x.unwrap()).collect();
    println!("No. of Sequences {}", rcds.len());
    let nseq = rcds.len();

    let now = std::time::Instant::now();

    let mut vidx : usize = 0;
    let mut pairwise_dcs: Vec<f64> = vec![0.0; (nseq*(nseq-1))/2];
    for x in 0..nseq{
        let seqx = rcds[x].seq();
        for y in (x+1)..nseq {
            let seqy = rcds[y].seq();
            // let tlen = tvec.len();

            let acsk = adyar_pairwise_acsk(&rargs, seqx, seqy);
            //println!("{}", acsk);

            pairwise_dcs[vidx] = acsk;
            vidx += 1;
        }
    // }
    }

    if rargs.show_progress == true {
        println!("D");
    }
    println!("Runtime for computing ACSK {} ms", now.elapsed().as_millis());
    let mut titles: Vec<String> = Vec::with_capacity(nseq);
    for x in 0..nseq {
        titles.push(String::from(rcds[x].id()))
    }
    write_output(&rargs, &titles, &pairwise_dcs)?;
    Ok(())
}

fn main() -> std::io::Result<()> {
    fn is_integer(v: String) -> Result<(), String> {
        if v.parse::<usize>().is_ok() { return Ok(()); }
        Err(String::from("The value is not valid integer"))
    }
    fn is_writable(v: String) -> Result<(), String> {
        let file = File::create(&v);  
        if file.is_ok() { return Ok(()); }
        Err(String::from("Can't open the output file"))
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
                            .arg(Arg::with_name("out_file")
                               .short("o")
                               .long("out_file")
                               .value_name("OUT_FILE")
                               .required(true)
                               .validator(is_writable)
                               .help("No. of Mismatches")
                               .takes_value(true))
                            .arg(Arg::with_name("r")
                               .short("r")
                               .long("use_rmq")
                               .help("Indicator to user rmq or not")
                               .takes_value(false))
                            .arg(Arg::with_name("p")
                               .short("p")
                               .long("show_progress")
                               .help("Indicator to show progress")
                               .takes_value(false))
                          .get_matches();

    
    let rargs: RunArgs = RunArgs {
        file_name: String::from(matches.value_of("INPUT").unwrap()),
        k: matches.value_of("mismatches").unwrap().parse::<usize>().unwrap(),
        out_file: String::from(matches.value_of("out_file").unwrap()),
        use_rmq: matches.occurrences_of("r") > 0,
        show_progress: matches.occurrences_of("p") > 0,
    };
    println!("Using input args: INPUT = {}, K = {}, Output = {}, Use RMQ {} ", rargs.file_name,
                 rargs.k, rargs.out_file, rargs.use_rmq);
    let tnow = std::time::Instant::now();
    run_adyar(rargs)?;
    println!("Total run time w. IO {} ms", tnow.elapsed().as_millis());
    Ok(())
}
