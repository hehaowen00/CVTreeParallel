#![feature(const_for)]
#![feature(const_mut_refs)]
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

const CODE: [i8; 26] = [ 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 ];
const LEN: usize = 6;
const AA_NUMBER: i64 = 20;
const EPSILON: f64 = 1e-010;

const M2: i64 = {
    let mut i: i64 = 1;
    let mut j = 0;
    while j < LEN - 2 {
        i = i * AA_NUMBER;
        j += 1;
    }
    i
};
const M1: i64 = M2 * AA_NUMBER;
const M: i64 = M1 * AA_NUMBER;

fn encode(ch: char) -> i8 {
    CODE[(ch as usize) - ('A' as usize)]
}

#[derive(Debug, PartialEq)]
struct Bacteria {
    count: i64,
    tv: Vec<f64>,
    ti: Vec<i64>,
    one_l: Vec<i64>,
    vector: Vec<i64>,
    second: Vec<i64>,
    total_l: i64,
    total: i64,
    indexs: i64,
    complement: i64,
}

impl Bacteria {
    pub fn new() -> Self {
        Self {
            count: 0,
            tv: Vec::new(),
            ti: Vec::new(),
            one_l: vec![0; AA_NUMBER as usize],
            second: vec![0; M1 as usize],
            vector: vec![0; M as usize],
            total:  0,
            total_l: 0,
            indexs: 0,
            complement: 0,
        }
    }

    pub fn init_buffer(
        &mut self,
        buf: &[u8]
    )
    {
        self.complement += 1;
        self.indexs = 0;
        for i in 0..LEN -1 {
            let enc: i8 = encode(buf[i] as char);
            self.one_l[enc as usize] += 1;
            self.total_l += 1;
            self.indexs = (self.indexs) * AA_NUMBER as i64 + enc as i64;
        }
        self.second[self.indexs as usize] += 1;
    }

    pub fn cont_buffer(
        &mut self,
        ch: u8
    )
    {
        let enc: i8 = encode(ch as char);
        self.one_l[enc as usize] += 1;
        self.total_l += 1;
        let index = self.indexs * AA_NUMBER as i64 + enc as i64;
        self.vector[index as usize] += 1;
        self.total += 1;
        self.indexs = (self.indexs % M2 as i64) * AA_NUMBER  as i64 + enc as i64;
        self.second[self.indexs as usize] += 1;
    }

    pub fn from_file(path: &str) -> Self {
        let mut bac = Bacteria::new();
        let f = File::open(path).unwrap(); 

        let mut reader = BufReader::new(f);
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).unwrap();

        let mut idx = 0;
        let mut b = bytes[idx];

        loop {
            if b == b'>' {
                while b != b'\n' {
                    idx += 1;
                    b = bytes[idx];
                }
                idx += 1;
                bac.init_buffer(&bytes[idx..idx + LEN - 1]);
            } else if b != b'\n' {
                bac.cont_buffer(b);
            }

            idx += 1;
            if idx >= bytes.len() {
                break;
            }
            b = bytes[idx];
        }

        bac
    }

    pub fn stochastic_compute(&self, i: i64) -> f64 {
        let p1 = self.second[(i / AA_NUMBER as i64) as usize] as f64 / (self.total + self.complement) as f64;
        let p2 = self.one_l[(i % AA_NUMBER as i64) as usize] as f64 / (self.total_l as f64);
        let p3 = self.second[(i % M1 as i64) as usize] as f64 / ((self.total + self.complement) as f64);
        let p4 = self.one_l[(i / M1 as i64) as usize]  as f64/ (self.total_l as f64);
        self.total as f64 * (p1 *p2 + p3 * p4) / 2.0
    }

    pub fn compare(&self, other: &Self) -> f64 {
        let mut correlation: f64 = 0.0;
        let mut vector_len1: f64 = 0.0;
        let mut vector_len2: f64 = 0.0;
        println!("{:?}  {} {}", self.one_l, self.indexs, self.total);
        return 0.0

//         println!("M {}", M);
//         for i in 0..M {
//             let stochastic1 = self.stochastic_compute(i as i64);
//             let  t1 = {
//                 if stochastic1 > EPSILON {
//                     (self.vector[i as usize]  as f64 - stochastic1) / stochastic1
//                 } else {
//                     0.0
//                 }
//             };
//             vector_len1 += t1 * t1;
//             // println!("t1 {}", t1);

//             let stochastic2 = other.stochastic_compute(i as i64);
//             let t2 = {
//                 if stochastic2 > EPSILON {
//                     (other.vector[i as usize] as f64 - stochastic2) / stochastic2
//                 } else {
//                     0.0
//                 }
//             };
//             // println!("t2 {}");

//             vector_len2 += t2 * t2;
//             correlation = correlation + t1 * t2;
//         }

//         println!("{}, {:?}, {:?}", correlation, vector_len1, vector_len2);

//         correlation /  vector_len1.sqrt() * vector_len2.sqrt()
    }
}

fn read_input_file() -> (usize, Vec<String>) {
    let  f = File::open("list.txt").unwrap();
    let reader = BufReader::new(f);
    let mut xs = Vec::new();

    let mut lines = reader.lines();
    let count = lines.next().unwrap().unwrap();
    for line in lines {
        xs.push(format!("data/{}.faa", line.unwrap()));
    }

    let count = usize::from_str_radix(&count, 10);
    (count.unwrap(), xs)
}

fn compare_all(count: usize, files: Vec<String>) {
    // let mut b = Vec::with_capacity(count);
    // for i in 0..count {
    //     println!("load {} of {}", i, count);
    // }
    for i in 0..count {
        let b1 = Bacteria::from_file(&files[0]);
        for j in i + 1.. count {
            let b2 = Bacteria::from_file(&files[1]);
            let correlation = b1.compare(&b2);
            print!("{:.2} {:.2} -> ", i, j);
            println!("{:.10}", correlation);
            break;
        }
        break;
    }
}

fn main() {
    let (count, files) = read_input_file();
    compare_all(count, files);
}

// Precompute all the pairs to be compared
// Load all files needed asynchronously
// As files are loaded, they are converted into bacteria instances
