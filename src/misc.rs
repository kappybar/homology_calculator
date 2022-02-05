pub fn betti_sphere(n: u32) -> Vec<u32> {
    let n: usize = n as usize;
    let mut s = vec![0; n + 1];
    s[0] = 1;
    s[n] = 1;
    s
}
