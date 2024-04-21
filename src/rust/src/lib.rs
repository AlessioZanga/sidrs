use std::{collections::VecDeque, hash::BuildHasherDefault, usize};

use indexmap::IndexSet;
use ndarray::prelude::*;
use rayon::prelude::*;
use rustc_hash::FxHasher;

pub type FxIndexSet<T> = IndexSet<T, BuildHasherDefault<FxHasher>>;

// Macro to get the parents of a vertex i in a graph G.
macro_rules! Pa {
    ($g:expr, $i:expr) => {
        $g.column($i)
            .iter()
            .enumerate()
            .filter_map(|(j, &a)| if a { Some(j) } else { None })
    };
}

// Macro to get the children of a vertex i in a graph G.
macro_rules! Ch {
    ($g:expr, $i:expr) => {
        $g.row($i)
            .iter()
            .enumerate()
            .filter_map(|(j, &a)| if a { Some(j) } else { None })
    };
}

// Macro to get the ancestors of a vertex i in a graph G.
macro_rules! An {
    ($g:expr, $i:expr) => {{
        let mut an = Vec::new();
        let mut to_be_visited = VecDeque::new();
        to_be_visited.push_back($i);
        while let Some(w) = to_be_visited.pop_front() {
            an.push(w);
            to_be_visited.extend(Pa!($g, w));
        }
        an.into_iter()
    }};
}

// Compute the path matrix of a DAG G.
fn _path(g: &Array2<bool>) -> Array2<bool> {
    // Get the number of vertices in G.
    let n = g.shape()[1];
    // Initialize the path matrix.
    let mut p = g.to_owned();

    // NOTE: We set the diagonal to true, as the path from a vertex to itself is always true.
    p.diag_mut().fill(true);

    // Compute the path matrix using the Floyd-Warshall algorithm.
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                p[[i, j]] = p[[i, j]] || (p[[i, k]] && p[[k, j]]);
            }
        }
    }

    p
}

// Get the reachable vertices from x given z in a graph G.
fn _reach(g: &Array2<bool>, x: usize, z: &[usize]) -> FxIndexSet<usize> {
    // Phase I - Get all ancestors of Z.
    let an_z: FxIndexSet<_> = z.iter().flat_map(|&z| An!(g, z)).collect();

    // Phase II - Traverse the active trail from X to Y.

    // Initialize the set of to-be-visited vertices.
    let mut to_be_visited = VecDeque::with_capacity(2 * g.shape()[1]);
    to_be_visited.push_back((x, true));
    // Initialize the set of visited vertices.
    let mut visited = FxIndexSet::<(usize, bool)>::default();
    visited.reserve(2 * g.shape()[1]);

    // While there are vertices to be visited.
    while let Some((w, d)) = to_be_visited.pop_front() {
        // Check if current vertex has already been visited.
        if visited.contains(&(w, d)) {
            continue;
        }
        // Add current vertex to visited set.
        visited.insert((w, d));
        // Trail up through W if W not in Z.
        if d && !z.contains(&w) {
            // Add parents of W to to-be-visited set.
            to_be_visited.extend(Pa!(g, w).map(|w| (w, true)));
            // Add children of W to to-be-visited set.
            to_be_visited.extend(Ch!(g, w).map(|w| (w, false)));
        // Trail down through W.
        } else if !d {
            // If W is not in Z, add children of W to to-be-visited set.
            if !z.contains(&w) {
                to_be_visited.extend(Ch!(g, w).map(|w| (w, false)));
            }
            // If W is in the ancestral set of Z, add parents of W to to-be-visited set.
            if an_z.contains(&w) {
                to_be_visited.extend(Pa!(g, w).map(|w| (w, true)));
            }
        }
    }

    // Return the set of visited vertices.
    visited.into_iter().map(|(w, _)| w).collect()
}

fn _sid_count(
    i: usize,
    j: usize,
    g_pa_i: &[usize],
    h_pa_i: &[usize],
    g: &Array2<bool>,
    p: &Array2<bool>,
    r: &FxIndexSet<usize>,
) -> usize {
    let g_ij_null = !p[[i, j]];
    let h_ij_null = h_pa_i.binary_search(&j).is_ok();

    if g_ij_null && h_ij_null {
        return 0;
    }

    if !g_ij_null && h_ij_null {
        return 1;
    }

    if g_pa_i == h_pa_i {
        return 0;
    }

    if g.row(i)
        .iter()
        .zip(p.column(j))
        .enumerate()
        .filter_map(|(k, (&a, &b))| if a && b { Some(k) } else { None })
        .any(|k| h_pa_i.iter().any(|&l| p[[k, l]]))
    {
        return 1;
    }

    if r.contains(&j) {
        return 1;
    }

    0
}

fn _sid(g: &Array2<bool>, h: &Array2<bool>) -> usize {
    // Initialize the structural intervention distance.
    let mut sid = 0;
    // Get the number of vertices in G.
    let n = g.shape()[1];
    // Compute the path matrix of G.
    let p = _path(g);
    // For each vertex i ...
    for i in 0..n {
        // Compute the parents of i in G and H.
        let g_pa_i: Vec<_> = Pa!(g, i).collect();
        let h_pa_i: Vec<_> = Pa!(h, i).collect();
        // Compute the vertices reachable from i in G.
        let r = _reach(g, i, &h_pa_i);
        // For each vertex j ...
        for j in (0..n).filter(|&j| j != i) {
            // Increment the SID count.
            sid += _sid_count(i, j, &g_pa_i, &h_pa_i, g, &p, &r);
        }
    }

    sid
}

fn _par_sid(g: &Array2<bool>, h: &Array2<bool>) -> usize {
    // Get the number of vertices in G.
    let n = g.shape()[1];
    // Compute the path matrix of G.
    let p = &_path(g);
    // For each vertex i ...
    (0..n)
        .into_par_iter()
        .flat_map(|i| {
            // Compute the parents of i in G and H.
            let g_pa_i: Vec<_> = Pa!(g, i).collect();
            let h_pa_i: Vec<_> = Pa!(h, i).collect();
            // Compute the vertices reachable from i in G.
            let r = _reach(g, i, &h_pa_i);
            // For each vertex j ...
            (0..n).into_par_iter().filter_map(move |j| {
                if j != i {
                    // Increment the SID count.
                    Some(_sid_count(i, j, &g_pa_i, &h_pa_i, g, p, &r))
                } else {
                    None
                }
            })
        })
        .sum()
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_path() {
        // Test with a 3x3 matrix.
        let g = array![
            [false, true, false],
            [false, false, true],
            [false, false, false]
        ];
        let p_true = array![
            [true, true, true],
            [false, true, true],
            [false, false, true]
        ];
        let p_pred = _path(&g);
        assert_eq!(p_true, p_pred);

        // Test with a 4x4 matrix.
        let g = array![
            [false, true, false, false],
            [false, false, true, false],
            [false, false, false, true],
            [false, false, false, false]
        ];
        let p_true = array![
            [true, true, true, true],
            [false, true, true, true],
            [false, false, true, true],
            [false, false, false, true]
        ];
        let p_pred = _path(&g);
        assert_eq!(p_true, p_pred);
    }

    #[test]
    fn test_d_sep() {
        let g = array![
            [false, true, false],
            [false, false, true],
            [false, false, false]
        ];
        assert!(_reach(&g, 0, &[1]).iter().eq(&[0, 1]));
        assert!(_reach(&g, 0, &[0]).iter().eq(&[0]));
        assert!(_reach(&g, 0, &[2]).iter().eq(&[0, 1, 2]));
        assert!(_reach(&g, 0, &[0, 1]).iter().eq(&[0]));
        assert!(_reach(&g, 0, &[0, 2]).iter().eq(&[0]));
        assert!(_reach(&g, 0, &[1, 2]).iter().eq(&[0, 1]));
    }

    #[test]
    #[rustfmt::skip]
    fn test_sid() {
        // Test from the paper example 1.
        let g = array![
            //  X1,    X2,    Y1,    Y2,    Y3
            [false,  true,  true,  true,  true], // X1
            [false, false,  true,  true,  true], // X2
            [false, false, false, false, false], // Y1
            [false, false, false, false, false], // Y2
            [false, false, false, false, false], // Y3
        ];
        let h_1 = array![
            //  X1,    X2,    Y1,    Y2,    Y3
            [false,  true,  true,  true,  true], // X1
            [false, false,  true,  true,  true], // X2
            [false, false, false,  true, false], // Y1
            [false, false, false, false, false], // Y2
            [false, false, false, false, false], // Y3
        ];
        let h_2 = array![
            //  X1,    X2,    Y1,    Y2,    Y3
            [false, false,  true,  true,  true], // X1
            [ true, false,  true,  true,  true], // X2
            [false, false, false, false, false], // Y1
            [false, false, false, false, false], // Y2
            [false, false, false, false, false], // Y3
        ];

        assert_eq!(_sid(&g, &h_1), 0);
        assert_eq!(_sid(&g, &h_2), 8);
    }

    #[test]
    #[rustfmt::skip]
    fn test_par_sid() {
        // Test from the paper example 1.
        let g = array![
            //  X1,    X2,    Y1,    Y2,    Y3
            [false,  true,  true,  true,  true], // X1
            [false, false,  true,  true,  true], // X2
            [false, false, false, false, false], // Y1
            [false, false, false, false, false], // Y2
            [false, false, false, false, false], // Y3
        ];
        let h_1 = array![
            //  X1,    X2,    Y1,    Y2,    Y3
            [false,  true,  true,  true,  true], // X1
            [false, false,  true,  true,  true], // X2
            [false, false, false,  true, false], // Y1
            [false, false, false, false, false], // Y2
            [false, false, false, false, false], // Y3
        ];
        let h_2 = array![
            //  X1,    X2,    Y1,    Y2,    Y3
            [false, false,  true,  true,  true], // X1
            [ true, false,  true,  true,  true], // X2
            [false, false, false, false, false], // Y1
            [false, false, false, false, false], // Y2
            [false, false, false, false, false], // Y3
        ];

        assert_eq!(_par_sid(&g, &h_1), 0);
        assert_eq!(_par_sid(&g, &h_2), 8);
    }
}

// Export the `sid` function to R.
use extendr_api::{prelude::*, ToVectorValue};

/// Convert an Robj to an Array2<bool>.
fn robj_to_array2_bool(obj: Robj) -> Array2<bool> {
    // Invoke `as.matrix`.
    let obj = R!("as.matrix({{obj}})").expect("Failed to invoke `as.matrix` on R object");
    // Get underlying shape and data vector.
    match obj.rtype() {
        // RMatrix<Rfloat> => Array2<bool>
        Rtype::Doubles => {
            let obj: RMatrix<Rfloat> = obj
                .as_matrix()
                .expect("Failed to convert to RMatrix<Rfloat>");
            Array::from_iter(obj.data().iter().map(|&x| x.to_real() > 0.)).into_shape(*obj.dim())
        }
        // RMatrix<Rint> => Array2<bool>
        Rtype::Integers => {
            let obj: RMatrix<Rint> = obj.as_matrix().expect("Failed to convert to RMatrix<Rint>");
            Array::from_iter(obj.data().iter().map(|&x| x.to_integer() > 0)).into_shape(*obj.dim())
        }
        // RMatrix<Rbool> => Array2<bool>
        Rtype::Logicals => {
            let obj: RMatrix<Rbool> = obj
                .as_matrix()
                .expect("Failed to convert to RMatrix<Rbool>");
            Array::from_iter(obj.data().iter().map(|&x| x.to_bool())).into_shape(*obj.dim())
        }
        _ => panic!(
            "Invalid R object class `{:?}` with type `{:?}`",
            obj.class(),
            obj.rtype()
        ),
    }
    .expect("Failed to convert to Array2<bool>")
}

/// Compute the SID between two adjacency matrices.
/// @export
#[extendr]
fn sid(g: Robj, h: Robj) -> usize {
    _sid(&robj_to_array2_bool(g), &robj_to_array2_bool(h))
}

/// Compute the SID between two adjacency matrices in parallel.
/// @export
#[extendr]
fn par_sid(g: Robj, h: Robj) -> usize {
    _par_sid(&robj_to_array2_bool(g), &robj_to_array2_bool(h))
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod sidrs;
    fn sid;
    fn par_sid;
}
