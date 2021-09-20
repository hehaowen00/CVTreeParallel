use std::io::Write;
use std::sync::Arc;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, AsyncReadExt, AsyncWriteExt, BufReader, stdout};
use tokio::sync::{RwLock, mpsc::Sender};

const AA_NUMBER: i64 = 20; // number of amino acids
const LEN: i64 = 6;
const EPSILON: f64 = 1e-010;
const CODE: [i8; 26] = [
    0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3,
];

const M2: i64 = AA_NUMBER.pow((LEN-2) as u32);
const M1: i64 = M2 * AA_NUMBER; // number of different possible 5-mers
const M: i64 = M1 * AA_NUMBER; // number of different possible 6-mers

const fn encode(ch: u8) -> i8 {
    CODE[(ch as u8 - 'A' as u8) as usize]
}

#[derive(Debug, Clone)]
struct Bacteria {
    one_l: Vec<i64>,
    indexs: i64,
    total: i64,
    total_l: i64,
    complement: i64,
    count: i64,
    cs: Vec<(i64, f64)>,
}

impl Bacteria {
    fn init_vectors() -> Self {
        Self {
            count: 0,
            one_l: vec![0; AA_NUMBER as usize],
            total: 0,
            total_l: 0,
            complement: 0,
            indexs: 0,
            cs: Vec::new(),
        }
    }

    fn init_buffer(&mut self, s: &[u8], second: &mut [f64]) {
        self.complement += 1;
        self.indexs = 0;
        for i in 0..LEN - 1 {
            let enc = encode(s[i as usize] as u8);
            self.one_l[enc as usize] += 1;
            self.total_l += 1;
            self.indexs = self.indexs * AA_NUMBER as i64 + enc as i64;
        }
        second[self.indexs as usize] += 1.0;
    }

    fn cont_buffer(&mut self, ch: u8, second: &mut [f64], vector: &mut [i64]) {
        let enc = encode(ch as u8);
        self.one_l[enc as usize] += 1;
        self.total_l += 1;
        let index = self.indexs * AA_NUMBER + enc as i64;
        vector[index as usize] += 1;
        self.total += 1;
        self.indexs = (self.indexs % M2) * AA_NUMBER + enc as i64;
        second[self.indexs as usize] += 1.0;
    }

    async unsafe fn read(&mut self, filename: &str, a: &mut [i64], b: &mut [f64]) {
        superluminal_perf::begin_event("file read start");
        let f = File::open(filename).await.unwrap();
        let mut s = BufReader::new(f);

        let vector = a;
        let (one_l_div_total, second) = b.split_at_mut(AA_NUMBER as usize);

        let mut buf: [u8; (LEN - 1) as usize] = [0; (LEN - 1) as usize];
        let mut ch: [u8; 1] = [0; 1];
        let mut sink = Vec::with_capacity(110);

        while let Ok(i) = s.read(&mut ch).await {
            if i == 0 {
                break;
            };
            if ch[0] == b'>' {
                s.read_until(b'\n', &mut sink).await;
                s.read_exact(&mut buf).await;
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

        for i in 0..AA_NUMBER {
            one_l_div_total[i as usize] = self.one_l[i as usize] as f64 / self.total_l as f64;
        }

        for i in 0..M1 {
            second[i as usize] = second[i as usize] / total_complement as f64;
        }

        self.count = 0;

        for i in 0..M {
            let i_div_m1 = i as i64 / M1 as i64;
            let i_mod_m1 = i as i64 - (M1 as i64 * i_div_m1);
            let i_div_aa_number = i / AA_NUMBER;
            let i_mod_aa_number = i - (AA_NUMBER * i_div_aa_number);

            let p1 = second[i_div_aa_number as usize];
            let p2 = one_l_div_total[i_mod_aa_number as usize];
            let p3 = second[i_mod_m1 as usize];
            let p4 = one_l_div_total[i_div_m1 as usize];
            let sto = (p1 * p2 + p3 * p4) * total_div_2;

            if sto > EPSILON {
                let r = (vector[i as usize] as f64 - sto) / sto;
                self.cs.push((i, r));
                self.count += 1;
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
            let (n1, t1) = self.cs[p1 as usize];
            let (n2, t2) = rhs.cs[p2 as usize];

            if n1 < n2 {
                v_len1 += t1 * t1;
                p1 += 1;
            } else if n2 < n1 {
                v_len2 += t2 * t2;
                p2 += 1;
            } else {
                p1 += 1;
                p2 += 1;

                v_len1 += t1 * t1;
                v_len2 += t2 * t2;
                correlation += t1 * t2;
            }
        }
        let total: f64 = self.cs[p1 as usize..].iter().map(|(_, v)| v * v).sum();
        v_len1 += total;
        
        let total: f64 = rhs.cs[p2 as usize..].iter().map(|(_, v)| v * v).sum();
        v_len2 += total;

        let res = correlation / (v_len1.sqrt() * v_len2.sqrt());
        superluminal_perf::end_event();
        res
    }
}

async fn read_input_file(tx: async_channel::Sender<(usize, String)>, fname: &str) {
    let f = File::open(fname).await.unwrap();
    let reader = BufReader::new(f);

    let mut lines = reader.lines();
    lines.next_line().await.unwrap().unwrap();
    let mut i = 0;
    while let Some(line) = lines.next_line().await.unwrap() {
        tx.send((i, format!("data/{}.faa", line))).await;
        i += 1;
    }
}

async fn load_bacteria(
    rx: async_channel::Receiver<(usize, String)>,
    notifier: Sender<(usize, String)>,
    bacterias: Arc<Vec<RwLock<Bacteria>>>
) {
    let mut a = vec![0; (M) as usize];
    let mut bx = vec![0.0; (M1 + AA_NUMBER) as usize];

    while let Ok((index, filename)) = rx.recv().await {
        let mut b = bacterias[index].write().await;
        unsafe {
            b.read(&filename, &mut a, &mut bx).await;
        }
        drop(b);

        notifier.send((index, filename)).await;

        a.clear();
        a.resize((M) as usize, 0);
        bx.clear();
        bx.resize((M1 + AA_NUMBER) as usize, 0.0);
    }
}


async fn compare(
    i: usize,
    j: usize,
    bacterias: Arc<Vec<RwLock<Bacteria>>>,
    out: Sender<(usize, usize, f64)>) {
    superluminal_perf::begin_event("compare one");
    let b1 = bacterias[i].read().await;
    let b2 = bacterias[j].read().await;
    let res = b1.compare(&b2);
    superluminal_perf::end_event(); 
    out.send((i, j, res)).await;
}

#[tokio::main]
async fn main() {
    let mut bacterias_ = Vec::with_capacity(41);
    for _ in 0..41 {
        bacterias_.push(RwLock::new(Bacteria::init_vectors()));
    }
    let bacterias = Arc::new(bacterias_);
    let (tx, rx) = async_channel::bounded(41);
    let (tx1, mut rx1) = tokio::sync::mpsc::channel(41);
    let (tx2, mut rx2) = tokio::sync::mpsc::channel(820);

    let mut done = Vec::with_capacity(41);
    let mut buf = Vec::with_capacity(27000);
    
    let num = 10; // number of cores + 2
    for _ in 0..num {
        tokio::spawn(load_bacteria(rx.clone(), tx1.clone(), bacterias.clone()));
    }
    drop(tx1);
    
    let t1 = std::time::Instant::now();
    superluminal_perf::begin_event("read input file");
    read_input_file(tx, "list.txt").await;
    superluminal_perf::end_event();

    tokio::spawn(async move {
        let mut buf = Vec::with_capacity(1100);
        while let Some((index, filename)) = rx1.recv().await {
            for j in &done {
                tokio::spawn(compare(index, *j, bacterias.clone(), tx2.clone()));
            }
            done.push(index);
            write!(&mut buf, "loaded {}, {}\n", index, filename);
        }
        stdout().write_all(&buf).await;
        drop(tx2);
    });

    while let Some((i, j, res)) = rx2.recv().await {
        if i < j {
            write!(&mut buf, "{:03} {:03} -> {:.10}\n", i, j, res);
        } else {
            write!(&mut buf, "{:03} {:03} -> {:.10}\n", j, i, res);
        }
    }
    std::io::stdout().write_all(&buf);
    let t2 = std::time::Instant::now();

    let diff = t2 - t1;
    println!("\n{} milliseconds\n", diff.as_millis());
}
