use homology_calculator::homology::{betti, Homology};

fn main() {
    let sphere1 = Homology::Sphere(1);
    let sphere2 = Homology::Sphere(1);

    let torus = Homology::Product(Box::new(sphere1), Box::new(sphere2));
    let betti_torus = betti(&torus);
    println!("{:?}", betti_torus);
}
