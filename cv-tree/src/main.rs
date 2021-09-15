#![feature(const_for)]
#![feature(const_mut_refs)]

use chashmap::CHashMap;
use std::sync::atomic::AtomicPtr;
use std::sync::Arc;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, AsyncReadExt, BufReader};

const AA_NUMBER: i64 = 20; // number of amino acids
const LEN: i64 = 6;
const EPSILON: f64 = 1e-010;
const CODE: [i8; 26] = [
    0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3,
];

const M2: i64 = {
    let mut i = 1;
    let mut j = 0;

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

#[derive(Debug, Clone)]
struct Bacteria {
    one_l: Vec<i64>, // [i64; AA_NUMBER as usize],
    indexs: i64,
    total: i64,
    total_l: i64,
    complement: i64,
    count: i64,
    tv: Vec<f64>,
    ti: Vec<i64>,
}

impl Bacteria {
    fn init_vectors() -> Self {
        Self {
            count: 0,
            one_l: [0; (AA_NUMBER as usize)].to_vec(),
            total: 0,
            total_l: 0,
            complement: 0,
            indexs: 0,
            tv: vec![],
            ti: vec![],
        }
    }

    fn init_buffer(&mut self, s: &[u8], second: &mut [i64]) {
        self.complement += 1;
        self.indexs = 0;
        for i in 0..LEN - 1 {
            let enc = encode(s[i as usize] as u8);
            self.one_l[enc as usize] += 1;
            self.total_l += 1;
            self.indexs = self.indexs * AA_NUMBER as i64 + enc as i64;
        }
        second[self.indexs as usize] += 1;
    }

    fn cont_buffer(&mut self, ch: u8, second: &mut [i64], vector: &mut [i64]) {
        let enc = encode(ch as u8);
        self.one_l[enc as usize] += 1;
        self.total_l += 1;
        let index = self.indexs * AA_NUMBER + enc as i64;
        vector[index as usize] += 1;
        self.total += 1;
        self.indexs = (self.indexs % M2) * AA_NUMBER + enc as i64;
        second[self.indexs as usize] += 1;
    }

    async fn read(&mut self, filename: &str, a: &mut [i64], b: &mut [f64]) {
        superluminal_perf::begin_event("file read start");
        let f = File::open(filename).await.unwrap();
        let mut s = BufReader::new(f);

        let (second, vector) = a.split_at_mut(M1 as usize);
        let (one_l_div_total, m2) = b.split_at_mut(AA_NUMBER as usize);
        let (second_div_total, t) = m2.split_at_mut(M1 as usize);

        let mut buf: [u8; (LEN - 1) as usize] = [0; (LEN - 1) as usize];
        let mut ch: [u8; 1] = [0; 1];
        let mut sink = Vec::with_capacity(110);

        while let Ok(i) = s.read(&mut ch).await {
            if i == 0 {
                break;
            };
            if ch[0] == b'>' {
                let _ = s.read_until(b'\n', &mut sink).await;
                // while ch[0] != b'\n' {
                //     s.read(&mut ch).await;
                // }
                let _ = s.read_exact(&mut buf).await;
                self.init_buffer(&buf, second);
            } else if ch[0] != b'\n' && ch[0] != b'\r' {
                self.cont_buffer(ch[0], second, vector);
            }

            sink.clear();
        }

        superluminal_perf::end_event();

        superluminal_perf::begin_event("precompute");
        let total_complement = self.total + self.complement;
        let total_div_2: f64 = self.total as f64 * 0.5;
        let mut i_mod_aa_number: i32 = 0;
        let mut i_div_aa_number: i32 = 0;
        let mut i_mod_m1: i64 = 0;
        let mut i_div_m1: i64 = 0;

        // let mut temp = [0.0; 20];

        // #[cfg(target_feature = "avx")]
        // unsafe {
        //     use std::arch::x86_64::*;
        //     for i in 0..AA_NUMBER / 4 {
        //         let start = i as usize * 4;
        //         let a = [
        //             self.one_l[start] as f64,
        //             self.one_l[start + 1] as f64,
        //             self.one_l[start + 2] as f64,
        //             self.one_l[start + 3] as f64
        //         ];
        //         let a = _mm256_loadu_pd(a.as_ptr());
        //         let b: [f64; 4] = std::mem::transmute([self.total_l as f64; 4]);
        //         let b = _mm256_loadu_pd(b.as_ptr());
        //         let c = _mm256_div_pd(a, b);
        //         let xs: [f64; 4] = std::mem::transmute(c);
        //         one_l_div_total[start] = xs[0];
        //         one_l_div_total[start + 1] = xs[1];
        //         one_l_div_total[start + 2] = xs[2];
        //         one_l_div_total[start + 3] = xs[3];
        //     }
        // }

        for i in 0..AA_NUMBER {
            one_l_div_total[i as usize] = self.one_l[i as usize] as f64 / self.total_l as f64;
        }

        for i in 0..M1 {
            second_div_total[i as usize] = second[i as usize] as f64 / total_complement as f64;
        }

        self.count = 0;

        for i in 0..M {
            let p1 = second_div_total[i_div_aa_number as usize];
            let p2 = one_l_div_total[i_mod_aa_number as usize];
            let p3 = second_div_total[i_mod_m1 as usize];
            let p4 = one_l_div_total[i_div_m1 as usize];
            let sto = (p1 * p2 + p3 * p4) * total_div_2;

            if i_mod_aa_number == AA_NUMBER as i32 - 1 {
                i_mod_aa_number = 0;
                i_div_aa_number += 1;
            } else {
                i_mod_aa_number += 1;
            }

            if i_mod_m1 == M1 - 1 {
                i_mod_m1 = 0;
                i_div_m1 += 1;
            } else {
                i_mod_m1 += 1;
            }

            if sto > EPSILON {
                t[i as usize] = (vector[i as usize] as f64 - sto) / sto;
                self.count += 1;
            } else {
                t[i as usize] = 0.0;
            }
        }

        self.ti.resize(self.count as usize, 0);
        self.tv.resize(self.count as usize, 0.0);

        let mut pos = 0;
        for i in 0..M {
            if t[i as usize] != 0.0 { 
                self.tv[pos] = t[i as usize];
                self.ti[pos] = i;
                pos += 1;
            }
        }
        superluminal_perf::end_event();
    }

    fn compare(&self, rhs: &Self) -> f64 {
        superluminal_perf::begin_event("compare start");
        let mut correlation: f64 = 0.0;
        let mut v_len1: f64 = 0.0;
        let mut v_len2: f64 = 0.0;

        let mut p1: i64 = 0;
        let mut p2: i64 = 0;

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
        }

        while p1 < self.count {
            // let n1 = self.ti[p1 as usize]; // does nothing
            let t1 = self.tv[p1 as usize];
            p1 += 1;
            v_len1 += t1 * t1;
        }

        while p2 < rhs.count {
            // let n2 = rhs.ti[p2 as usize]; // does nothing
            let t2 = rhs.tv[p2 as usize];
            p2 += 1;
            v_len2 += t2 * t2;
        }

        let res = correlation / (v_len1.sqrt() * v_len2.sqrt());
        superluminal_perf::end_event();
        res
    }
}

async fn read_input_file(tx: async_channel::Sender<(usize, String)>, fname: &str) {
    let f = File::open(fname).await.unwrap();
    let reader = BufReader::new(f);

    let mut lines = reader.lines();
    let _ = lines.next_line().await.unwrap().unwrap();
    let mut i = 0;
    while let Some(line) = lines.next_line().await.unwrap() {
        let _ = tx.send((i, format!("data/{}.faa", line))).await;
        i += 1;
    }
}

async fn load_bacteria(
    rx: async_channel::Receiver<(usize, String)>,
    notifier: async_channel::Sender<usize>,
    bacterias: Arc<CHashMap<usize, Bacteria>>
) {
    let mut a = vec![0; (M1 + M) as usize];
    let mut bx = vec![0.0; (M1 + M + AA_NUMBER) as usize];

    while let Ok((index, filename)) = rx.recv().await {
        let mut b = bacterias.remove(&index).unwrap();
        b.read(&filename, &mut a, &mut bx).await;
        print!("loaded {}, {}\n", index, filename);
        
        bacterias.insert(index, b);
        notifier.send(index).await;

        a.clear();
        a.resize((M1 + M) as usize, 0);
        bx.clear();
        bx.resize((M1 + M + AA_NUMBER) as usize, 0.0);
    }
}

async fn compare(i: usize, j: usize, bacterias: Arc<CHashMap<usize, Bacteria>>) {
    superluminal_perf::begin_event("compare one");
    let b1 = bacterias.get(&i).unwrap();
    let b2 = bacterias.get(&j).unwrap();
    superluminal_perf::end_event(); 
    print!("{:03} {:03} -> {}\n", i, j, b1.compare(&b2));
}

#[tokio::main]
async fn main() {
    let mut handles = Vec::with_capacity(8);
    let mut bacterias = Arc::new(CHashMap::with_capacity(64));
    let (tx, rx) = async_channel::bounded(41);
    let (tx1, rx1) = async_channel::bounded(41);

    for i in 0..41 {
        bacterias.insert(i, Bacteria::init_vectors());
    }

    let num = 10;

    for i in 0..num  {
        tokio::spawn(load_bacteria(rx.clone(), tx1.clone(), bacterias.clone()));
    }

    let t1 = std::time::Instant::now();

    superluminal_perf::begin_event("read input file sse");
    let names = read_input_file(tx, "list.txt").await;
    superluminal_perf::end_event();

    // load_bacteria(rx, tx1, bacterias.clone()).await;

    // while let Some(h) = handles.pop() {
    //     h.await;
    // }

    let mut done = Vec::with_capacity(41);
    while let Ok(index) = rx1.recv().await {
        for j in &done {
            let h = tokio::spawn(compare(index, *j, bacterias.clone()));
            handles.push(h);
        }
        done.push(index);
        if done.len() == 41 {
            break;
        }
    }

    // for i in 0..41 {
    //     for j in i + 1..41 {
    //         let h = tokio::spawn(compare(i, j, bacterias.clone()));
    //         handles.push(h);
    //     }
    // }

    while let Some(h) = handles.pop() {
        h.await;
    }

    let t2 = std::time::Instant::now();

    let diff = t2 - t1;
    println!("{} milliseconds", diff.as_millis());
}
