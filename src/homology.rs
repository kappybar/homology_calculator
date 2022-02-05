use crate::misc;
use crate::simplical_complex;
use std::cmp::max;

pub enum Homology {
    SimplicalComplex(Vec<Vec<u32>>),
    Sphere(u32),
    Union(Box<Homology>, Box<Homology>),
    Product(Box<Homology>, Box<Homology>),
}

pub fn betti(homology: &Homology) -> Vec<u32> {
    match homology {
        Homology::SimplicalComplex(sc) => simplical_complex::betti_simplical_complex(sc),
        Homology::Sphere(n) => misc::betti_sphere(*n),
        Homology::Union(h1, h2) => betti_union(h1, h2),
        Homology::Product(h1, h2) => betti_product(h1, h2),
    }
}

fn betti_union(homology1: &Box<Homology>, homology2: &Box<Homology>) -> Vec<u32> {
    let h1 = betti(&homology1);
    let h2 = betti(&homology2);
    let dim = max(h1.len(), h2.len());
    let mut h = vec![0; dim];
    for (i, h1i) in h1.iter().enumerate() {
        h[i] += h1i;
    }
    for (i, h2i) in h2.iter().enumerate() {
        h[i] += h2i;
    }
    h
}

fn betti_product(homology1: &Box<Homology>, homology2: &Box<Homology>) -> Vec<u32> {
    let h1 = betti(&homology1);
    let h2 = betti(&homology2);
    let dim = h1.len() + h2.len() - 1;
    let mut h = vec![0; dim];
    for (i, h1i) in h1.iter().enumerate() {
        for (j, h2j) in h2.iter().enumerate() {
            h[i + j] += h1i * h2j;
        }
    }
    h
}
