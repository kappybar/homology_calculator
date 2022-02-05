use homology_calculator::homology::{betti, Homology};

fn main() {
    let example = Homology::SimplicalComplex(vec![vec![1, 2, 3], vec![2, 3, 4], vec![1, 4]]);
    let betti_example = betti(&example);
    println!("{:?}", betti_example);
}
