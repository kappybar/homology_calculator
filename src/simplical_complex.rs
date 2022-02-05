use itertools::Itertools;
use std::cmp::{max, min};
use std::mem;

pub fn betti_simplical_complex(simplical_complex: &Vec<Vec<u32>>) -> Vec<u32> {
    let max_dim = {
        let mut dim = 0;
        for simlex in simplical_complex {
            dim = max(dim, simlex.len());
        }
        dim
    };
    let simpleces = extract_simleces(simplical_complex, max_dim);
    let c: Vec<usize> = simpleces.iter().map(|x| x.len()).collect();
    let mut h: Vec<usize> = vec![0; max_dim];

    for k in 0..(max_dim - 1) {
        // calculate h_k
        let mut representation_matrix = vec![vec![0; c[k + 1]]; c[k]];
        for (i, simplex) in simpleces[k + 1].iter().enumerate() {
            for j in 0..simplex.len() {
                let mut boundary = Vec::new();
                for l in 0..simplex.len() {
                    if l != j {
                        boundary.push(simplex[l]);
                    }
                }

                let index = simpleces[k]
                    .iter()
                    .enumerate()
                    .find(|(_i, x)| **x == boundary)
                    .unwrap()
                    .0;
                representation_matrix[index][i] = if j % 2 == 0 { 1 } else { -1 };
            }
        }
        h[k] = smith_normal_form(&mut representation_matrix);
    }

    let mut betti: Vec<u32> = vec![0; max_dim];
    for i in 0..max_dim {
        if i == 0 {
            betti[i] = (c[i] - h[i]) as u32;
        } else {
            betti[i] = (c[i] - h[i] - h[i - 1]) as u32;
        }
    }
    while let Some(last) = betti.last() {
        if *last != 0 {
            break;
        } else {
            betti.pop();
        }
    }

    betti
}

fn extract_simleces(simplical_complex: &Vec<Vec<u32>>, max_dim: usize) -> Vec<Vec<Vec<u32>>> {
    let mut simpleces: Vec<Vec<Vec<u32>>> = vec![Vec::new(); max_dim];
    for simplex in simplical_complex.iter() {
        for bit in 1..(1 << simplex.len()) {
            let subset: Vec<u32> = (0..simplex.len())
                .filter(|&i| (bit >> i & 1) == 1)
                .map(|i| simplex[i])
                .collect();
            let count = subset.len() - 1;
            simpleces[count].push(subset);
        }
    }
    for ndim_simpleces in &mut simpleces {
        for i in 0..ndim_simpleces.len() {
            ndim_simpleces[i].sort();
        }
    }
    simpleces
        .into_iter()
        .map(|x| x.into_iter().unique().collect())
        .collect()
}

// find absolute minimum element (not equal 0) in [k,n] * [k,n] submatrix
fn find_abs_minimum_element(matrix: &Vec<Vec<i32>>, k: usize) -> (usize, usize) {
    let mut index_i = k;
    let mut index_j = k;
    let mut min_elem = i32::MAX;
    for i in k..matrix.len() {
        for j in k..matrix[i].len() {
            if matrix[i][j].abs() < min_elem && matrix[i][j] != 0 {
                index_i = i;
                index_j = j;
                min_elem = matrix[i][j].abs();
            }
        }
    }
    (index_i, index_j)
}

fn swap_row(matrix: &mut Vec<Vec<i32>>, i1: usize, i2: usize) {
    let mut tmp = Vec::new();
    mem::swap(&mut matrix[i1], &mut tmp);
    mem::swap(&mut matrix[i2], &mut tmp);
    mem::swap(&mut matrix[i1], &mut tmp);
}

fn swap_column(matrix: &mut Vec<Vec<i32>>, j1: usize, j2: usize) {
    for row in matrix.iter_mut() {
        row.swap(j1, j2);
    }
}

fn column_transform(matrix: &mut Vec<Vec<i32>>, i: usize, k: usize) {
    let r = matrix[i][k] / matrix[k][k];
    for j in k..matrix[i].len() {
        matrix[i][j] -= matrix[k][j] * r;
    }
}

fn row_transform(matrix: &mut Vec<Vec<i32>>, j: usize, k: usize) {
    let r = matrix[k][j] / matrix[k][k];
    for i in k..matrix.len() {
        matrix[i][j] -= matrix[i][k] * r;
    }
}

fn smith_normal_form(matrix: &mut Vec<Vec<i32>>) -> usize {
    let size = min(matrix.len(), matrix[0].len());
    for k in 0..size {
        loop {
            loop {
                let (index_i, index_j) = find_abs_minimum_element(matrix, k);

                if index_i == k && index_j == k {
                    break;
                }

                // swap row
                if index_i != k {
                    swap_row(matrix, index_i, k);
                }

                // swap column
                if index_j != k {
                    swap_column(matrix, index_j, k);
                }

                // k-column
                // minimize nondividable element
                for i in (k + 1)..matrix.len() {
                    if matrix[i][k] % matrix[k][k] != 0 {
                        column_transform(matrix, i, k);
                    }
                }

                // k-row
                // minimize nondividable element
                for j in (k + 1)..matrix[k].len() {
                    if matrix[k][j] % matrix[k][k] != 0 {
                        row_transform(matrix, j, k);
                    }
                }
            }

            if matrix[k][k] == 0 {
                return k;
            }

            // k-column
            for i in (k + 1)..matrix.len() {
                column_transform(matrix, i, k);
            }

            // k-row
            for j in (k + 1)..matrix[k].len() {
                row_transform(matrix, j, k);
            }

            let mut flag = true;
            for i in (k + 1)..matrix.len() {
                for j in (k + 1)..matrix[i].len() {
                    if matrix[i][j] % matrix[k][k] != 0 {
                        flag = false;
                        matrix[i][k] = matrix[k][k];
                        let r = matrix[i][j] / matrix[i][k];
                        for ii in k..matrix.len() {
                            matrix[ii][j] -= r * matrix[ii][k];
                        }
                        break;
                    }
                }
                if !flag {
                    break;
                }
            }

            if flag {
                matrix[k][k] = matrix[k][k].abs();
                break;
            }
        }
    }
    size
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        // smith_normal_form
        let mut matrix_test = vec![vec![2, 4, 4], vec![-6, 6, 12], vec![10, 4, 16]];
        let rank = smith_normal_form(&mut matrix_test);
        let matrix_correct = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 156]];
        assert_eq!(rank, 3);
        assert_eq!(matrix_test, matrix_correct);
        // betti_simplical_complex
        let simplical_complex = vec![vec![1, 2, 3], vec![2, 3, 4], vec![1, 4]];
        let betti_test = betti_simplical_complex(&simplical_complex);
        let betti_correct = vec![1, 1];
        assert_eq!(betti_test, betti_correct);
    }
}
