use homology_calculator::homology::{betti, Homology};

fn main() {
    let sphere = Homology::Sphere(1);
    let regular_tetrahedron = Homology::SimplicalComplex(vec![
        vec![1, 2, 3],
        vec![2, 3, 4],
        vec![1, 3, 4],
        vec![1, 2, 4],
    ]);

    let union = Homology::Union(Box::new(sphere), Box::new(regular_tetrahedron));
    let betti_union = betti(&union);
    println!("{:?}", betti_union);
}
