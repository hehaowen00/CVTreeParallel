#![feature(const_for)]
#![feature(const_mut_refs)]
use std::{fs::File, io::{BufReader, Read}};

const AA_NUMBER: i64 = 20;
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
};
const M1: i64 = M2 * AA_NUMBER;
const M: i64 = M1 * AA_NUMBER;

const fn encode(ch: u8) -> i8 {
    CODE[ch as usize - 'A' as usize]
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
}

impl Bacteria {
    fn init_vectors() -> Self {
        Self {
            vector: [0; (M as usize)].to_vec(),
            second: [0; (M1 as usize)].to_vec(),
            one_l: [0; (AA_NUMBER as usize)].to_vec(),
            total: 0,
            total_l: 0,
            complement: 0,
            indexs: 0,
        }
    }

    fn init_buffer(&mut self, s: &[u8]) {
        self.complement += 1;
        self.indexs = 0;
        for i in 0..LEN - 1 {
            let enc = encode(s[i as usize]);
            self.one_l[enc as usize] += 1;
            self.total_l += 1;
            self.indexs = self.indexs * AA_NUMBER as i64 + enc as i64;
        }
        self.second[self.indexs as usize] += 1;
    }

    fn cont_buffer(&mut self, ch: u8) {
        let enc = encode(ch);
        self.one_l[enc as usize] += 1;
        self.total_l += 1;
        let index = self.indexs * AA_NUMBER + enc as i64;
        self.vector[index as usize] += 1;
        self.total += 1;
        self.indexs = (self.indexs % M2) * AA_NUMBER + enc as i64;
        self.second[self.indexs as usize] += 1;
    }

    fn read(filename: &str) -> Self {
        let f = File::open(filename).unwrap();
        let mut reader = BufReader::new(f);
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).unwrap();

        let mut bac = Self::init_vectors();
        let mut idx = 0;
        let mut temp = buf[idx];

        loop {
            if temp == b'>' {
                while temp != b'\n' {
                    idx += 1;
                    temp = buf[idx];
                }
                idx += 1;
                print!("{:?}", std::str::from_utf8(&buf[idx .. idx + (LEN - 1) as usize]));
                bac.init_buffer(&buf[idx..idx + (LEN - 1) as usize]);
                idx =  idx + LEN as usize - 1;
            } else if temp != b'\n' {
print!("{}", temp as char);
                bac.cont_buffer(temp);
            }

            if idx + 1 >= buf.len() {
                break;
            }
            idx += 1;
            temp = buf[idx];
        }
        bac
    }
}

fn main() {
    // unsafe {
    //     M2 = 1;
    //     for _ in 0..LEN - 2 {
    //         M2 *= AA_NUMBER;
    //     }
    //     M1 = M2 * AA_NUMBER;
    //     M = M1 * AA_NUMBER;
    // }
    // unsafe {
    //     println!("M2 {} M1 {} M {}", M2, M1, M);
    //     for c in 'A'..'Z' {
    //         println!("{}", encode(c as u8));
    //     }
    // }
    // println!("AA {} LEN {} EPSI {}", AA_NUMBER, LEN, EPSILON);

    let bacteria = Bacteria::read("../data/AcMNPV.faa");
    // println!("{:?}", bacteria);
}
