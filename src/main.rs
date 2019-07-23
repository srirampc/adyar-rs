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

pub struct RunArgs {
    pub k: usize,
    pub file_name: String,
}

pub struct GST<'s, T> where 
  T: std::marker::Copy + num::Integer + std::fmt::Debug {
    pub stx: &'s [u8],
    pub sty: &'s [u8],
    pub txt: &'s [u8],
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


// int verify_acsk_at(RunArgs& args, GST& gst,
// 				   int pos, int match_pos,
// 				   int max_length){
pub fn verify_acsk_at<T>(gst: &GST<T>, pos: T,
                         match_pos: T, max_length: T, k: usize) -> usize
where 
  T: std::marker::Copy + 
     num::Integer + num::FromPrimitive + num::ToPrimitive +
     std::fmt::Debug {
    // int	count = 0;
    let mut count: usize = 0;
    let mut upos: usize = T::to_usize(&pos).unwrap();
    let mut umatch_pos: usize = T::to_usize(&match_pos).unwrap();
    let mut max_length: usize = T::to_usize(&max_length).unwrap();
	// // because lcpk = total length - k
	// max_length += args.k;
    max_length += k;
	// while(max_length){
    while max_length > 0 {
		// if(gst.text[pos] == gst.text[match_pos]){
        if gst.txt[upos] == gst.txt[umatch_pos]{
			// pos++; match_pos++;
            upos += 1; umatch_pos += 1;
		// }
		// else{
        } else {
			// count++;
			// pos++; match_pos++;
            count += 1; upos += 1; umatch_pos += 1;
		// }
        }
		// max_length--;
        max_length -= 1;
	// }
    }
	// if(count > args.k){
	// 	return count;
	// } else {
	// 	return 0;
	// }
    count
// }
}

// void verify_acsk(RunArgs& args, GST& gst,
// 				 std::vector<int>& lcpk_kmacs,
// 				 std::vector<int>& max_match){
pub fn verify_acsk<T>(gst: &GST<T>, lcpk: &LCPK<T>, k: usize) -> usize where 
  T: std::marker::Copy + 
     num::Integer + num::FromPrimitive + num::ToPrimitive +
     std::fmt::Debug {
	// // prints if # missmatches > k for all i in text 
	// int n = gst.text.length() + 1;
	// int num_err = 0;
    let n = gst.txt.len();
    let mut num_err: usize = 0;
	// for(int i=2; i<n; i++){ // skip 0 and 1 in SA
    for i in 2..n {
		// int count = verify_acsk_at(args, gst,
		// 						   gst.SA[i],
		// 						   gst.SA[max_match[i]],
		// 						   lcpk_kmacs[i]);
        let max_match = std::cmp::max(lcpk.left_match[i], lcpk.right_match[i]).
                            to_usize().unwrap();
        let count = verify_acsk_at(gst, gst.gsa.sarray[i], gst.gsa.sarray[max_match],
                                   lcpk.match_length[i], k);
		// if(count > 0){
		// 	std::cout << i << ": " << count << std::endl;
		// 	num_err++;
		// }
        if count != k {
            println!("{} : {}", i, count);
            num_err += 1;
        }
	// }
    }
	// std::cout << "Number of errors :" << num_err << std::endl;
    println!("No. of errors {}", num_err);
    num_err
// }
}

// void find_matches(RunArgs& args, GST& gst,
//                   std::vector<int>& match_length,
// 				  std::vector<int>& left_match,
// 				  std::vector<int>& right_match){
pub fn find_matches<T>(gst: &GST<T>) -> LCPK<T> where 
  T: std::marker::Copy + num::Integer + num::FromPrimitive + std::fmt::Debug {
    
    // int n = gst.text.length() + 1;
	// match_length.resize(n, 0);
	// left_match.resize(n, 0);
	// right_match.resize(n, 0);
    let n = gst.txt.len();
    let mut match_length: Vec<T> = vec![T::zero(); n];
    let mut left_match: Vec<T> = vec![T::zero(); n];
    let mut right_match: Vec<T> = vec![T::zero(); n];
    let lcpk_length: Vec<T> = vec![T::zero(); n];
	// int m = args.y.length();
    let m = gst.sty.len();
    let eff_len = T::from_usize(n - m - 2).unwrap();
	// int min=0, last_match=0;
    let mut min: T = T::zero(); let mut last_match: T = T::zero();
	// for(int i=1; i < n-1; i++) {
    for i in 2..(n-1) {
		// if ((gst.SA[i] > n-m-2 && gst.SA[i+1] > n-m-2) || (gst.SA[i] < n-m-2 && gst.SA[i+1] < n-m-2)) {
        if (gst.gsa.sarray[i] > eff_len && gst.gsa.sarray[i+1] > eff_len) || 
            (gst.gsa.sarray[i] < eff_len && gst.gsa.sarray[i+1] < eff_len) {
			// if (gst.LCP[i+1] <= min)
			// 	min = gst.LCP[i+1];
            if gst.glcp[i+1] <= min {
                min = gst.glcp[i+1];
            }
			// match_length[i+1] = min;
            match_length[i+1] = min;
			// left_match[i+1] = last_match;
            left_match[i+1] = last_match;
		// }
		// else {
        } else {
			// min = gst.LCP[i+1];
            min = gst.glcp[i+1];
			// match_length[i+1] = gst.LCP[i+1];
            match_length[i+1] = gst.glcp[i+1];
			// left_match[i+1] = i; //stores the position in second string where the match occures
			// last_match = i;
            last_match = T::from_usize(i).unwrap();
            left_match[i+1] = last_match;
    	//}
        }
	// }
    }

	// min=0; last_match=0;
    min = T::zero(); last_match = T::zero();
	// for(int i=n-1; i >= 1; i--) {
    for i in (n-1)..0 {
		// if ((gst.SA[i] > n-m-2 && gst.SA[i-1] > n-m-2) || (gst.SA[i] < n-m-2 && gst.SA[i-1] < n-m-2)) {
        if (gst.gsa.sarray[i] > eff_len && gst.gsa.sarray[i-1] > eff_len) ||
            (gst.gsa.sarray[i] < eff_len && gst.gsa.sarray[i-1] < eff_len) {
			// if (gst.LCP[i] <= min)
			// 	min = gst.LCP[i];
            if gst.glcp[i] <= min {
                min = gst.glcp[i];
            }
			// if(min>match_length[i-1]){
				// match_length[i-1] = min;
				// right_match[i-1] = last_match;
				// left_match[i-1] = 0;
			// }
			// else if(min<match_length[i-1])
			// 	right_match[i-1]=0;
			// else
			// 	right_match[i-1]=last_match;
            if min > match_length[i-1] {
				match_length[i-1] = min;
				right_match[i-1] = last_match;
				left_match[i-1] = T::zero();
            } else if min<match_length[i-1] {
				right_match[i-1]= T::zero();
            } else {
				right_match[i-1]=last_match;
            }
		// }
		// else {
        } else {
			// min = gst.LCP[i];
			// last_match = i;
			min = gst.glcp[i];
			last_match = T::from_usize(i).unwrap();
			// if(min>match_length[i-1]){
				// match_length[i-1] = min;
				// right_match[i-1] = last_match;
				// left_match[i-1] = 0;
			// }
			// else if(min<match_length[i-1])
			// 	right_match[i-1] = 0;
			// else
			// 	right_match[i-1] = last_match;
			if min > match_length[i-1] {
				match_length[i-1] = min;
				right_match[i-1] = last_match;
				left_match[i-1] = T::zero();
			}
			else if min < match_length[i-1] {
				right_match[i-1] = T::zero();
			} else {
				right_match[i-1] = last_match;
            }
		// }
        }
	// }
    }

    LCPK::<T>{
        match_length: match_length,
        left_match: left_match,
        right_match: right_match,
        lcpk_length: lcpk_length
    }
// }
}

// int forward_match(std::string const& T, int pos1, int pos2, int k, int m, int n){
pub fn forward_match(txt: &[u8], pos1: usize, pos2: usize, 
                     k: usize, m: usize, n: usize) -> usize {
	// int lcp = 0;
    let mut lcp: usize = 0;
    let mut rk : i32 = k as i32;
    let mut p1 = pos1;
    let mut p2 = pos2;
    println!("n {} m {}", n, m);
	// while(k>=0 && pos1< n-m-2 && pos2 < n){
    while rk >= 0 && p1 < (n - m - 2) && p2 < n {
		// if(T[pos1] == T[pos2]){
		// 	lcp++; pos1++; pos2++;
		// }
		// else {
		// 	pos1++; pos2++; k--;
		// }
        if txt[p1] == txt[p2] {
            lcp += 1; p1 += 1; p2 += 1;
        } else{
            p1 += 1; p2 += 1; rk -= 1;
        }
	// }
    }
	// return lcp;
    lcp
// }
}

#[derive(Debug, PartialEq)]
pub enum MatchDirectionFlag{
    MatchNone,
    MatchLeft,
    MatchRight
}

pub fn match_with_second<T>(args: & RunArgs, gst: & GST<T>,
                        lcpk: &LCPK<T>, pos1: usize, i: usize) -> (MatchDirectionFlag, usize, usize)
    where
    T: std::marker::Copy + 
         num::Integer + num::FromPrimitive + num::ToPrimitive +
        std::fmt::Debug {
    let n = gst.txt.len();
    let m = gst.sty.len();
    let mut pos2: usize = 0;
    let mut maxl: usize = 0;
    let mut flag: MatchDirectionFlag = MatchDirectionFlag::MatchNone;
    let mut p: usize; let mut lcp: usize;
    if lcpk.left_match[i] > T::zero() && pos1 < n-m-2 {
        p = lcpk.left_match[i].to_usize().unwrap();
        while gst.glcp[p+1] >= lcpk.match_length[i] && p > 0 {
            if gst.gsa.sarray[p].to_usize().unwrap() > n-m-2 {
                lcp = forward_match(gst.txt,
                                    pos1-1,
                                    (gst.gsa.sarray[p] +
                                        lcpk.match_length[i]).to_usize().unwrap(),
                                    args.k, n-m-2, n);
                if lcp > maxl {
                    maxl = lcp;
                    pos2 = p;
                    flag = MatchDirectionFlag::MatchLeft;
                }
            }
            p -= 1;
        }
    }
    if lcpk.right_match[i] > T::zero() && pos1 < n-m-2 {
        p = lcpk.right_match[i].to_usize().unwrap();
        while gst.glcp[p] >= lcpk.match_length[i] && p < n {
            if gst.gsa.sarray[p].to_usize().unwrap() > n-m-2 {
                lcp = forward_match(gst.txt,
                                    pos1-1,
                                    (gst.gsa.sarray[p] + 
                                        lcpk.match_length[i]).to_usize().unwrap(),
                                    args.k, n-m-2, n);
                if lcp > maxl {
                    maxl =lcp;
                    pos2 = p;
                    flag = MatchDirectionFlag::MatchRight;
                }
            }
            p += 1;
        }
    }
    (flag, pos2, maxl)
}

pub fn match_with_first<T>(args: & RunArgs, gst: & GST<T>,
                        lcpk: &LCPK<T>, pos1: usize, i: usize) -> (MatchDirectionFlag, usize, usize)
    where
    T: std::marker::Copy + 
         num::Integer + num::FromPrimitive + num::ToPrimitive +
        std::fmt::Debug {

    let n = gst.txt.len();
    let m = gst.sty.len();

    let mut pos2: usize = n;
    let mut maxl: usize = 0;
    let mut flag: MatchDirectionFlag = MatchDirectionFlag::MatchNone;
    let mut p: usize; let mut lcp: usize;

    if lcpk.left_match[i] > T::zero() && pos1 < n {
        p = lcpk.left_match[i].to_usize().unwrap();
        while gst.glcp[p+1] >= lcpk.match_length[i] && p > 0 {
            if gst.gsa.sarray[p].to_usize().unwrap() < n-m-2 {
                lcp = forward_match(gst.txt,
                                    (gst.gsa.sarray[p]+
                                        lcpk.match_length[i]).to_usize().unwrap(), 
                                    pos1-1,
                                    args.k, n-m-2, n);
                if lcp > maxl {
                    maxl = lcp;
                    pos2 = p;
                    flag = MatchDirectionFlag::MatchLeft;
                }
            }
            p -= 1;
        }
    }
    if lcpk.right_match[i] > T::zero() && pos1 < n {
        p = lcpk.right_match[i].to_usize().unwrap();
        while gst.glcp[p] >= lcpk.match_length[i] && p < n {
            if gst.gsa.sarray[p].to_usize().unwrap() < n-m-2 {
                lcp = forward_match(gst.txt,
                                    (gst.gsa.sarray[p] +
                                        lcpk.match_length[i]).to_usize().unwrap(),
                                    pos1-1,
                                    args.k, n-m-2, n);
                if lcp > maxl {
                    maxl = lcp;
                    pos2 = p;
                    flag = MatchDirectionFlag::MatchRight;
                }
            }
            p += 1;
        }
    }
    (flag, pos2, maxl)
}

pub fn compute_lcpk_kmacs<T>(args: & RunArgs, gst: & GST<T>, lcpk: &mut LCPK<T>) -> i32
where 
  T: std::marker::Copy + 
     num::Integer + num::FromPrimitive + num::ToPrimitive +
     std::fmt::Debug {
    let status: i32 = 0;
    let n = gst.txt.len();
    let m = gst.sty.len();
    assert!(lcpk.lcpk_length.len() == n);

    for i in 0..n {
        lcpk.lcpk_length[i] = lcpk.match_length[i];
        let pos1 = gst.gsa.sarray[i].to_usize().unwrap()  +
                    lcpk.match_length[i].to_usize().unwrap();

		// forward matching char by char
		let (flag, pos2, maxl) = if gst.gsa.sarray[i].to_usize().unwrap() < n-m-2 {
            match_with_second(args, gst, lcpk, pos1, i)
		} else if gst.gsa.sarray[i].to_usize().unwrap() > n-m-2 {
            match_with_first(args, gst, lcpk, pos1, i)
		} else {
            (MatchDirectionFlag::MatchNone, 0, 0)
        };
    
        if flag == MatchDirectionFlag::MatchLeft {
            lcpk.left_match[i] = T::from_usize(pos2).unwrap();
            lcpk.right_match[i] = T::zero();
        }
        else if flag == MatchDirectionFlag::MatchRight {
            lcpk.left_match[i] = T::zero();
            lcpk.right_match[i] = T::from_usize(pos2).unwrap();
        }
        lcpk.lcpk_length[i] = lcpk.match_length[i] + T::from_usize(maxl).unwrap();

    }

    status
}


// void compute_acsk(RunArgs& args, GST& gst,
// 				  std::vector<int>& match_length,
// 				  std::vector<int>& left_match,
// 				  std::vector<int>& right_match){
pub fn compute_acsk<T>(args: & RunArgs, gst: & GST<T>, lcpk: &mut LCPK<T>) -> f64
    where 
    T: std::marker::Copy + 
        num::Integer + num::FromPrimitive + num::ToPrimitive +
        std::fmt::Debug {

	// lcpk based on k-macs's heuristic
	// std::vector<int> lcpk_kmacs;

	compute_lcpk_kmacs(args, gst, lcpk);

	// to calculate and print acsk based in k-macs
	//compute_acsk_kmacs(args, gst , lcpk_kmacs);

    let n = gst.txt.len();
	let m = gst.sty.len();
	let mut max_match: Vec<T> = vec![T::zero(); n];

	for i in 0..n {
		if lcpk.left_match[i] < lcpk.right_match[i] {
			max_match[i] = lcpk.right_match[i];
		} else {
			max_match[i] = lcpk.left_match[i];
        }
	}
	// to verify k-macs algorithm
	//verify_acsk(args, gst, lcpk_kmacs, max_match);

	// backward and forward search and storage of missmatches for each i
    let mut mismatch_position: Vec<usize> = vec![n; 2*args.k+2];
	for i in 2..n { // skip i = 0 and i = 1 as first two chars are string terminatros
		let mut kidx = args.k + 1; 
		// initialize pos1 , pos2
		let mut pos1: usize = gst.gsa.sarray[i].to_usize().unwrap(); 
        let mut pos2: usize = gst.gsa.sarray[max_match[i].to_usize().unwrap()].to_usize().unwrap();
		//backward matching

		while ((pos1 < n-m-2 && pos2 > n-m-2) || (pos1 > n-m-2 && pos2 < n-m-2)) && kidx > 0 {
			if !(gst.txt[pos1] == gst.txt[pos2]) {
				if (pos1 < n-m-2 && pos2 < n-m-2) || (pos1 > n-m-2 && pos2 > n-m-2) { //checks for '$' crossover
					break;
				}
                mismatch_position[kidx-1] = pos1 + 1;
				kidx -= 1;
			}
            if pos1 == 0 || pos2 == 0 {
                break;
            }
            pos1 -= 1; pos2 -= 1;
		}

		//forward matching
		pos1 = gst.gsa.sarray[i].to_usize().unwrap(); 
        pos2 = gst.gsa.sarray[max_match[i].to_usize().unwrap()].to_usize().unwrap();
		kidx = args.k + 1;

		while ((pos1 < n-m-2 && pos2 > n-m-2) || (pos1 > n-m-2 && pos2 < n-m-2)) && kidx > 0 && pos1 < n && pos2 < n {
			if !(gst.txt[pos1] == gst.txt[pos2]) {
				if (pos1 < n-m-2 && pos2 < n-m-2) || (pos1 > n-m-2 && pos2 > n-m-2) {//checks for '$' crossover
					break;
				}
                mismatch_position[2*args.k - (kidx-1) + 1] = pos1-1; 
				kidx -= 1;
			}
            pos1 += 1; pos2 += 1;
		}
		// for each missmatch compair sk[i] and [k+j+1]th - [j]th
		// elements of mismatch_position array replace if less
    }

	let mut max_match_new: Vec<T> = max_match.clone(); // vec![T::zero(); n];
	// for i in 0..n { max_match_new[i] = max_match[i]; }
    for i in 2..n { // skip i = 0 and i = 1 as first two chars are string terminatros
		for j in 0..args.k {
			if mismatch_position[j] == n || mismatch_position[j+args.k+1] == n {
				continue;
			}
            let mut substring_len = mismatch_position[args.k+j+1] - mismatch_position[j]+1;
			//subtract k from substring to comply to acs definition
			substring_len -= args.k;
            let mkidx = gst.gisa[mismatch_position[j]].to_usize().unwrap();
			if substring_len > lcpk.lcpk_length[mkidx].to_usize().unwrap() {
				lcpk.lcpk_length[mkidx] = T::from_usize(substring_len).unwrap();
				// modify max_match_new[] to reflect new lcpk_kmacs value,
				// needed for verify_acsk() function 
				let origin_position = gst.gsa.sarray[max_match[i].to_usize().unwrap()].to_usize().unwrap();
				let offset = gst.gsa.sarray[i].to_usize().unwrap() - mismatch_position[j];
				max_match_new[mkidx] = gst.gisa[origin_position - offset];
			}
		}
	}
	// to verify phase 1 of adyar algorithm
	//verify_acsk(args, gst, lcpk_kmacs, max_match_new);
	// modify sk-array so that sk[i] = x implies sk[i+1] >= x-1
	for i in 1..n {
        let skidx = gst.gisa[gst.gsa.sarray[i].to_usize().unwrap() + 1].to_usize().unwrap();
		if lcpk.lcpk_length[i] > lcpk.lcpk_length[skidx] {
			lcpk.lcpk_length[skidx] = lcpk.lcpk_length[i] - T::one();
			// modify max_match_new for phase 2 to account for positions
			// that were not caught by mismatch_postion array.
			max_match_new[skidx] = 
                gst.gisa[gst.gsa.sarray[max_match_new[i].to_usize().unwrap()].to_usize().unwrap()+1];
		}
	}
	// to verify adyar algorithm
	//verify_acsk(args, gst, lcpk_kmacs, max_match_new);
	let mut s1: Vec<usize> = vec![0; n-m-2];
    let mut s2: Vec<usize> = vec![0; m];
	for i in 1..(n-1) {
        let saidx = gst.gsa.sarray[i].to_usize().unwrap();
		if saidx >= n-m-1 {
			s2[saidx-(n-m-1)] = lcpk.lcpk_length[i].to_usize().unwrap();
		}
        if saidx <= n-m-3 {
			s1[saidx] = lcpk.lcpk_length[i].to_usize().unwrap();
	    }
    }

	// calculate avgs1 and svgs2 and d_acs
	let mut avg_s1: f64=0.0; let mut avg_s2: f64=0.0;
	for i in 0..(n-m-2) {avg_s1 += s1[i] as f64;}
	avg_s1 = avg_s1/(n-m-2) as f64;
	for i in 0..m {avg_s2 += s2[i] as f64;}
	avg_s2 = avg_s2/m as f64;

	let d_acs: f64 =
     ((( (n-m-2) as f64).log10()/(2.0*avg_s2)) + ((m as f64).log10()/(2.0*avg_s1))) - 
        ((((n-m-2) as f64).log10()/(((n-m-2)  as f64))) + ((m as f64).log10()/(m as f64)));
	//std::cout << "ACS for k= " << args.k << " : " << d_acs << std::endl;
	// args.d_vector.push_back(d_acs);
    d_acs
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
            tvec.push('@' as u8);
            tvec.append(&mut seqy.to_vec());
            tvec.push('$' as u8);
            // let tlen = tvec.len();

            let gsa = sufsort::SA::<i32>::new(&tvec);
            let gisa = sufsort::construct_isa(&gsa.sarray);
            let glcp = lcp::construct_lcp_from_sa(gsa.txt, &gsa.sarray, &gisa);
            let grmq = rmq::RMQ::<i32, usize>::new(&glcp);

            let gst: GST<i32> = GST::<i32> {
                stx: &seqx, sty: &seqy, txt: &tvec, 
                gsa: &gsa, gisa: &gisa, glcp: &glcp, grmq: &grmq
            };

            let mut lcpk = find_matches(&gst);
            // println!("Matches {:?}", fmn.match_length);
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
