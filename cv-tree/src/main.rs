#![feature(const_for)]
#![feature(const_mut_refs)]
use std::{fs::{File}, io::{Read, BufReader, BufRead}, process::exit};

const AA_NUMBER: i64 = 20; // number of amino acids
const LEN: i64 = 6;
const EPSILON: f64 = 1e-010;
const CODE: [i8; 26] = [ 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 ];

const M2: i64 = {
    let mut i = 1;
    let mut j =0;

    while j < LEN - 2 {
        i *= AA_NUMBER;
        j += 1;
    }
    i
}; // number of different possible 4-mers
const M1: i64 = M2 * AA_NUMBER; // number of different possible 5-mers
const M: i64 = M1 * AA_NUMBER; // number of different possible 6-mers

const fn encode(ch: u8) -> i8 {
    CODE[(ch as u8 - 'A' as u8) as usize]
}

#[derive(Debug)]
struct Bacteria {
    second: Vec<i64>,
    one_l: Vec<i64>,// [i64; AA_NUMBER as usize],
    indexs: i64,
    total: i64,
    total_l: i64,
    complement: i64,
    vector: Vec<i64>,
    count: i64,
    tv: Vec<f64>,
    ti: Vec<i64>,
}

impl Bacteria {
    fn init_vectors() -> Self {
        Self {
            count: 0,
            vector: [0; (M as usize)].to_vec(),
            second: [0; (M1 as usize)].to_vec(),
            one_l: [0; (AA_NUMBER as usize)].to_vec(),
            total: 0,
            total_l: 0,
            complement: 0,
            indexs: 0,
            tv: vec![],
            ti: vec![],
        }
    }

    fn init_buffer(&mut self, s: &[u8]) {
        self.complement += 1;
        self.indexs = 0;
        for i in 0..LEN - 1 {
            let enc = encode(s[i as usize] as u8);
            self.one_l[enc as usize] += 1;
            self.total_l += 1;
            self.indexs = self.indexs * AA_NUMBER as i64 + enc as i64;
        }
        self.second[self.indexs as usize] += 1;
    }

    fn cont_buffer(&mut self, ch: u8) {
        let enc = encode(ch as u8);
        self.one_l[enc as usize] += 1;
        self.total_l += 1;
        let index = self.indexs * AA_NUMBER + enc as i64;
        self.vector[index as usize] += 1;
        self.total += 1;
        self.indexs = (self.indexs % M2) * AA_NUMBER + enc as i64;
        self.second[self.indexs as usize] += 1;
    }

    fn read(&mut self, filename: &str) {
        let mut f = File::open(filename).unwrap();
        let mut chars = String::new();
        f.read_to_string(&mut chars).unwrap();

        let mut buf: [u8; (LEN-1) as usize] = [0; (LEN-1) as usize];
        let s: Vec<char> = chars.chars().collect();
        let mut s = s.iter();

        while  let Some(mut ch) = s.next() {
            if *ch == '>' {
                while *ch != '\n' {
                    ch = s.next().unwrap();
                }
                for i in 0..(LEN -1) as usize {
                    buf[i] = *s.next().unwrap() as u8;
                }
                self.init_buffer(&buf);
            } else if ch.is_whitespace() == false {
                self.cont_buffer(*ch as u8);
            }
        }

        let total_complement = self.total + self.complement;
        let total_div_2: f64 = self.total as f64 * 0.5;
        let mut i_mod_aa_number: i32 = 0;
        let mut i_div_aa_number: i32 = 0;
        let mut i_mod_m1: i64 = 0;
        let mut i_div_m1: i64 = 0;

        let mut one_l_div_total = [0.0; AA_NUMBER as usize].to_vec();
        for i in 0..AA_NUMBER {
            one_l_div_total[i as usize] = self.one_l[i as usize] as f64 / self.total_l as f64;
        }
        
        let mut second_div_total = [0.0; M1 as usize].to_vec();
        for i in 0..M1 {
            second_div_total[i as usize] = self.second[i as usize] as f64 / total_complement as f64;
        }

        self.count = 0;
        let mut t = [0.0; M as usize].to_vec();

        for i in 0..M {
            let p1 = second_div_total[i_div_aa_number as usize]; 
            let p2 = one_l_div_total[i_mod_aa_number as usize];
            let p3 = second_div_total[i_mod_m1 as usize];
            let p4 = one_l_div_total[i_div_m1 as usize];
            let sto  = (p1 * p2 + p3 * p4) * total_div_2;

            if i_mod_aa_number == AA_NUMBER as i32 - 1 {
                i_mod_aa_number = 0;
                i_div_aa_number += 1;
            } else {
                i_mod_aa_number += 1;
            }

            if i_mod_m1 == M1-1 {
                i_mod_m1 = 0;
                i_div_m1 += 1;
            } else {
                i_mod_m1 += 1;
            }

            if sto > EPSILON {
                t[i as usize] = (self.vector[i as usize] as f64 - sto) / sto;
                self.count += 1;
            } else {
                t[i as usize] = 0.0;
            }
        }

        drop(second_div_total);
        self.vector = vec![];
        self.second = vec![]; 

        for _ in 0..self.count {
            self.tv.push(0.0);
            self.ti.push(0);
        }

        let mut pos = 0;
        for i in 0..M {
            if t[i as usize] != 0.0 {
                self.tv[pos] = t[i as usize];
                self.ti[pos] = i;
                pos += 1;
            }
        }
    }

    fn compare(&self, rhs: &Self) -> f64 {
        let mut correlation: f64 = 0.0;
        let mut v_len1: f64 = 0.0;
        let mut v_len2: f64 = 0.0;

        let mut p1: i64 = 0;
        let mut p2: i64 = 0;
        let mut i = 0;

        while p1 < self.count && p2 < rhs.count {
            let n1 = self.ti[p1 as usize];
            let n2 = rhs.ti[p2 as usize];

            if n1 < n2 {
                let t1 = self.tv[p1 as usize];
                v_len1 += t1 * t1;
                p1 += 1;
            } else if n2 < n1 {
                let t2 = rhs.tv[p2 as usize];
                v_len2 += t2 * t2; 
                p2 += 1;
            } else {
                let t1 = self.tv[p1 as usize];
                p1 += 1;

                let t2 = rhs.tv[p2 as usize];
                p2 += 1;

                v_len1 += (t1 * t1);
                v_len2 += (t2 * t2);
                correlation += t1 * t2;
            }

            i += 1;
        }

        while p1 < self.count {
            let n1 = self.ti[p1 as usize];
            let t1 = self.tv[p1 as usize];
            p1 += 1;
            v_len1 += t1 * t1;
        }

        while p2 < rhs.count {
            let n2 = rhs.ti[p2 as usize];
            let t2 = rhs.tv[p2 as usize];
            p2 += 1;
            v_len2 += t2 * t2;
        }

        correlation / (v_len1.sqrt() * v_len2.sqrt())
    }

//     fn stochastic(&self, i: i64) -> f64 {
//         let total_complement: f64 = (self.total + self.complement) as f64;
//         let p1: f64 = self.second[(i / AA_NUMBER) as usize] as f64 / (total_complement);
//         let p2: f64 = self.one_l[(i % AA_NUMBER) as usize] as f64 / (self.total_l as f64);
//         let p3: f64 = self.second[(i % M1) as usize] as f64 / total_complement;
//         let p4: f64 = self.one_l[(i / M1) as  usize] as f64 / (self.total_l as f64);
//         return self.total as f64 * (p1 * p2 + p3 * p4) / 2.0;
//     }

//     fn compare(&self, other: &Self) -> f64 {
//         let mut correlation: f64 = 0.0;
//         let mut v_len1: f64 = 0.0;
//         let mut v_len2: f64 = 0.0;

//         for i in 0..M {
//             let sto1 = self.stochastic(i as i64);
//             let t1 = if sto1 > EPSILON {
//                 (self.vector[i as usize] as f64 - sto1) / sto1
//             } else {
//                 0.0
//             };

//             v_len1 += (t1 * t1);

//             let sto2 = other.stochastic(i as i64);
//             let t2 = if sto2 > EPSILON {
//                 (other.vector[i as usize] as f64 - sto2) / sto2
//             } else {
//                 0.0
//             };

//             v_len2 += (t2 * t2);

//             correlation = correlation + t1 * t2;
//             if correlation == f64::NAN {
//                 correlation = 0.0;
//             }

//            if v_len1 == f64::NAN {
//                v_len1 = 0.0;
//            }

//            if v_len2 == f64::NAN {
//                v_len2 = 0.0;
//            }
            
//         }
//         correlation / (v_len1.sqrt() * v_len2.sqrt())
//     }
}

fn compare_all(names: Vec<String>) {
    let mut bacterias = Vec::with_capacity(names.len());
    for i in 0..names.len() {
        println!("load {} of {}", i + 1, names.len());
        let mut b = Bacteria::init_vectors();
        b.read(&names[i]);
        bacterias.push(b);
    }
}

fn read_input_file(fname: &str) -> Vec<String> {
    let  f = File::open(fname).unwrap();
    let reader = BufReader::new(f);
    let mut xs = Vec::new();

    let mut lines = reader.lines();
    let count = lines.next().unwrap().unwrap();
    for line in lines {
        xs.push(format!("data/{}.faa", line.unwrap()));
    }

    xs
}

fn main() {
//     let mut b1 = Bacteria::init_vectors();
//     b1.read("../data/AcMNPV.faa");
//     let mut b2 = Bacteria::init_vectors();
//     b2.read("../data/AdhoNPV.faa");

//     println!("{}", b1.compare(&b2));
    let names = read_input_file("list.txt");
    compare_all(names);
}
