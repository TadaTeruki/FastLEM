use std::collections::BinaryHeap;
use terrain_graph::edge_attributed_undirected::EdgeAttributedUndirectedGraph;

use crate::core::{
    traits::Site,
    units::{Elevation, Length},
};

/// Tree structure for representing the flow of water.
///  - `next` is the next site of each site in the flow.
pub struct StreamTree {
    pub next: Vec<usize>,
}

struct RidgeElement {
    index: usize,
    dist: Length,
}

impl RidgeElement {
    fn evaluate(&self) -> f64 {
        self.dist
    }
}

impl PartialEq for RidgeElement {
    fn eq(&self, other: &Self) -> bool {
        self.evaluate() == other.evaluate()
    }
}

impl Eq for RidgeElement {}

impl Ord for RidgeElement {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.evaluate().partial_cmp(&self.evaluate()).unwrap()
    }
}

impl PartialOrd for RidgeElement {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl StreamTree {
    /// Constructs a stream tree from a given terrain data.
    pub fn construct<S: Site>(
        sites: &[S],
        elevations: &[Elevation],
        graph: &EdgeAttributedUndirectedGraph<Length>,
        outlets: &[usize],
    ) -> Self {
        let num = sites.len();

        // `is_outlet` is a table that indicates whether a site is an outlet or not.
        let is_outlet = Self::create_outlet_table(sites, outlets);

        // `next` is the next site of each site in the flow.
        // at this point, the stream tree can create lakes: a root of a stream tree not connected to an outlet.
        let next = Self::construct_initial_stream_tree(num, elevations, graph, &is_outlet);

        // `subroot` is the root of each site in the flow. lakes are not removed yet.
        let (subroot, has_lake) = Self::find_roots_with_lakes(num, &is_outlet, &next);

        // if there are no lakes, stream tree is already complete
        if !has_lake {
            return StreamTree { next };
        }

        // remove lakes from the stream tree
        let next = Self::remove_lakes_from_stream_tree(&next, num, graph, outlets, &subroot);

        StreamTree { next }
    }

    fn create_outlet_table<S: Site>(sites: &[S], outlets: &[usize]) -> Vec<bool> {
        let mut is_outlet = vec![false; sites.len()];
        outlets.iter().for_each(|&i| {
            is_outlet[i] = true;
        });
        is_outlet
    }

    fn construct_initial_stream_tree(
        num: usize,
        elevations: &[Elevation],
        graph: &EdgeAttributedUndirectedGraph<Length>,
        is_outlet: &[bool],
    ) -> Vec<usize> {
        let mut next: Vec<usize> = (0..num).collect();

        (0..num).for_each(|i| {
            if is_outlet[i] {
                return;
            }

            let mut steepest_slope = 0.0;
            graph.neighbors_of(i).iter().for_each(|ja| {
                let j = ja.0;
                if elevations[i] > elevations[j] {
                    let distance = ja.1;
                    let down_hill_slope = (elevations[i] - elevations[j]) / distance;
                    if down_hill_slope > steepest_slope {
                        steepest_slope = down_hill_slope;
                        next[i] = j;
                    }
                }
            });
        });

        next
    }

    fn find_roots_with_lakes(num: usize, is_outlet: &[bool], next: &[usize]) -> (Vec<usize>, bool) {
        let mut subroot: Vec<Option<usize>> = (0..num)
            .map(|i| if is_outlet[i] { Some(i) } else { None })
            .collect();

        let mut has_lake = false;

        (0..num).for_each(|i| {
            if subroot[i].is_some() {
                return;
            }
            let mut iv = i;
            while subroot[iv].is_none() && iv != next[iv] {
                iv = next[iv];
            }

            let ir = {
                if subroot[iv].is_none() {
                    has_lake = true;
                    Some(iv)
                } else {
                    subroot[iv]
                }
            };

            let mut iv = i;
            while subroot[iv].is_none() && iv != next[iv] {
                subroot[iv] = ir;
                iv = next[iv];
            }
            subroot[iv] = ir;
        });

        (subroot.iter().map(|&r| r.unwrap()).collect(), has_lake)
    }

    fn remove_lakes_from_stream_tree(
        next: &[usize],
        num: usize,
        graph: &EdgeAttributedUndirectedGraph<Length>,
        outlets: &[usize],
        subroot: &[usize],
    ) -> Vec<usize> {
        // final roots of the stream tree
        let mut root: Vec<Option<usize>> = vec![None; num];
        let mut ridgestack: BinaryHeap<RidgeElement> = BinaryHeap::with_capacity(num);
        outlets.iter().for_each(|&outlet| {
            root[outlet] = Some(outlet);
            ridgestack.push(RidgeElement {
                index: outlet,
                dist: 0.0,
            });
        });

        // remove lakes
        let mut visited: Vec<bool> = vec![false; num];
        let mut next = next.to_owned();

        while let Some(element) = ridgestack.pop() {
            let i = element.index;

            if visited[i] {
                continue;
            }

            graph
                .neighbors_of(i)
                .iter()
                .enumerate()
                .for_each(|(_, ja)| {
                    let j = ja.0;
                    if visited[j] {
                        return;
                    }

                    if root[subroot[j]].is_none() {
                        let mut k = j;
                        let mut nk = i;
                        loop {
                            if next[k] != k {
                                // flip flow
                                let tmp = next[k];
                                next[k] = nk;
                                nk = k;
                                k = tmp;
                            } else {
                                break;
                            }
                        }
                        next[k] = nk;
                        root[subroot[j]] = root[subroot[i]];
                    }

                    let distance = ja.1;
                    ridgestack.push(RidgeElement {
                        index: j,
                        dist: distance,
                    });
                });
            root[i] = root[subroot[i]];
            visited[i] = true;
        }

        next
    }
}
