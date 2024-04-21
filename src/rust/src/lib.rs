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
fn _reach_non_directed(g: &Array2<bool>, x: usize, z: &FxIndexSet<usize>) -> FxIndexSet<usize> {
    // Remove outgoing edges from X.
    let mut g = g.to_owned();
    g.row_mut(x).fill(false);

    // Phase I - Get all ancestors of Z.
    let an_z: FxIndexSet<_> = z.iter().flat_map(|&z| An!(g, z)).collect();

    // Phase II - Traverse the active trail from X to Y.

    // Initialize the set of to-be-visited vertices.
    let mut to_be_visited = VecDeque::with_capacity(2 * g.shape()[1]);
    to_be_visited.push_back((x, true));
    // Initialize the set of visited vertices.
    let mut visited = FxIndexSet::<(usize, bool)>::default();
    visited.reserve(2 * g.shape()[1]);
    // Initialize the set of reachable vertices.
    let mut reachable = FxIndexSet::<usize>::default();
    reachable.reserve(g.shape()[1]);

    // While there are vertices to be visited.
    while let Some((w, d)) = to_be_visited.pop_front() {
        // Check if current vertex has already been visited.
        if visited.contains(&(w, d)) {
            continue;
        }
        // Add current vertex to visited set.
        visited.insert((w, d));
        // If W is not in Z, add W to reachable set.
        if !z.contains(&w) {
            reachable.insert(w);
        }
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
    reachable
}

fn _sid_count(
    x: usize,
    y: usize,
    g_ch_x: &[usize],
    h_pa_x: &FxIndexSet<usize>,
    p: &Array2<bool>,
    r: &FxIndexSet<usize>,
) -> usize {
    // Check if Y is a parent of X in H.
    if h_pa_x.contains(&y) {
        // Check if Y is a descendant of X in G.
        if p[[x, y]] {
            // Increment the SID count.
            return 1;
        }
    } else {
        // Check if Y is reachable from X in G from a non-directed path.
        if r.contains(&y) {
            // Increment the SID count.
            return 1;
        }
        // Get the vertices W onh_pa_x all directed paths from X to Y.
        let mut w = g_ch_x.iter().filter(|&&w| p[[w, y]]);
        // Check if there is a vertex W s.t. Z is a descendant of W.
        if w.any(|&w| h_pa_x.iter().any(|&z| w == z || p[[w, z]])) {
            // Increment the SID count.
            return 1;
        }
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
    for x in 0..n {
        // Compute the parents of i in G and H.
        let g_ch_x: Vec<_> = Ch!(g, x).collect();
        let h_pa_x: FxIndexSet<_> = Pa!(h, x).collect();
        // Compute the vertices reachable from i in G via non-directed paths.
        let r = _reach_non_directed(g, x, &h_pa_x);
        // For each vertex j ...
        for y in (0..n).filter(|&y| x != y) {
            // Increment the SID count.
            sid += _sid_count(x, y, &g_ch_x, &h_pa_x, &p, &r);
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
        .flat_map(|x| {
            // Compute the parents of i in G and H.
            let g_ch_x: Vec<_> = Ch!(g, x).collect();
            let h_pa_x: FxIndexSet<_> = Pa!(h, x).collect();
            // Compute the vertices reachable from i in G via non-directed paths.
            let r = _reach_non_directed(g, x, &h_pa_x);
            // For each vertex j ...
            (0..n).into_par_iter().filter_map(move |y| {
                if x != y {
                    // Increment the SID count.
                    Some(_sid_count(x, y, &g_ch_x, &h_pa_x, p, &r))
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
            [false, true, true],
            [false, false, true],
            [false, false, false]
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
            [false, true, true, true],
            [false, false, true, true],
            [false, false, false, true],
            [false, false, false, false]
        ];
        let p_pred = _path(&g);
        assert_eq!(p_true, p_pred);
    }

    #[test]
    fn test_reach_non_directed() {
        // X -> Y, Z -> X, Z -> Y.
        let g = array![
            // X, Y, Z
            [false, true, false],  // X
            [false, false, false], // Y
            [true, true, false],   // Z
        ];

        assert!(_reach_non_directed(&g, 0, &[].into_iter().collect())
            .iter()
            .eq(&[0, 2, 1]));
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

        assert_eq!(0, _sid(&g, &h_1));
        assert_eq!(8, _sid(&g, &h_2));
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

        assert_eq!(0, _par_sid(&g, &h_1));
        assert_eq!(8, _par_sid(&g, &h_2));
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
