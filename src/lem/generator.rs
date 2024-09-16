use rand::{rngs::StdRng, Rng, SeedableRng};
use std::marker::PhantomData;
use thiserror::Error;

use crate::{
    core::{
        parameters::TopographicalParameters,
        traits::{Model, Site},
        units::{Length, Step},
    },
    lem::drainage_basin::DrainageBasin,
    lem::stream_tree,
};

/// The default value of the exponent `m` for calculating stream power.
const DEFAULT_M_EXP: f64 = 0.5;

#[derive(Error, Debug)]
pub enum GenerationError {
    #[error("The number of topographical parameters must be equal to the number of sites")]
    InvalidNumberOfParameters,
    #[error("You must set topographical parameters before generating terrain")]
    ParametersNotSet,
    #[error("You must set `TerrainModel` before generating terrain")]
    ModelNotSet,
}

/// Provides methods for generating terrain.
///
/// ### Required properties
///  - `model` is the vector representation of the terrain network.
///  - `parameters` is the topographical parameters of sites. Each parameter contains the uplift rates, erodibilities, base elevations and maximum slopes (see [TopographicalParameters] for details).
/// ### Optional properties
///  - `max_iteration` is the maximum number of iterations. If not set, the iterations will be repeated until the elevations of all sites are stable.
///
#[derive(Clone)]
pub struct TerrainGenerator<S, M, T>
where
    S: Site,
    M: Model<S, T>,
{
    model: Option<M>,
    parameters: Option<Vec<TopographicalParameters>>,
    max_iteration: Option<Step>,
    _phantom: PhantomData<(S, T)>,
}

impl<S, M, T> Default for TerrainGenerator<S, M, T>
where
    S: Site,
    M: Model<S, T>,
{
    fn default() -> Self {
        Self {
            model: None,
            parameters: None,
            max_iteration: None,
            _phantom: PhantomData,
        }
    }
}

impl<S, M, T> TerrainGenerator<S, M, T>
where
    S: Site,
    M: Model<S, T>,
{
    /// Set the model that contains the set of sites.
    pub fn set_model(mut self, model: M) -> Self {
        self.model = Some(model);
        self
    }

    /// Set the topographical parameters of sites. See [TopographicalParameters] about the parameters.
    pub fn set_parameters(mut self, parameters: Vec<TopographicalParameters>) -> Self {
        self.parameters = Some(parameters);
        self
    }

    /// Set the maximum number of iterations.
    ///
    /// The iteration(loop) for calculating elevations will be stopped when the number of iterations reaches `max_iteration`.
    /// If not set, the iterations will be repeated until the elevations of all sites are stable.
    pub fn set_max_iteration(mut self, max_iteration: Step) -> Self {
        self.max_iteration = Some(max_iteration);
        self
    }

    /// Generate terrain.
    pub fn generate(self) -> Result<T, GenerationError> {
        let model = {
            if let Some(model) = &self.model {
                model
            } else {
                return Err(GenerationError::ModelNotSet);
            }
        };

        let (num, sites, areas, graph, default_outlets) = (
            model.num(),
            model.sites(),
            model.areas(),
            model.graph(),
            model.default_outlets(),
        );

        let parameters = {
            if let Some(parameters) = &self.parameters {
                if parameters.len() != num {
                    return Err(GenerationError::InvalidNumberOfParameters);
                }
                parameters
            } else {
                return Err(GenerationError::ParametersNotSet);
            }
        };

        let m_exp = DEFAULT_M_EXP;

        let outlets = {
            let outlets = parameters
                .iter()
                .enumerate()
                .filter(|(_, param)| param.is_outlet)
                .map(|(i, _)| i)
                .collect::<Vec<_>>();
            if outlets.is_empty() {
                default_outlets.to_vec()
            } else {
                outlets
            }
        };

        let mut rng: StdRng = SeedableRng::seed_from_u64(0);
        let mut elevations = parameters
            .iter()
            .map(|a| a.base_elevation + rng.gen::<f64>() * f64::EPSILON)
            .collect::<Vec<_>>();

            loop {

                let mut changed = false;
                if step < 1 {

                    let stream_tree =
                    stream_tree::StreamTree::construct(sites, &elevations, graph, &outlets);

                    let mut drainage_areas: Vec<f64> = areas.to_vec();
                    let mut response_times = vec![0.0; num];

                    // calculate elevations for each drainage basin
                    outlets.iter().for_each(|&outlet| {
                        // construct drainage basin
                        let drainage_basin = DrainageBasin::construct(outlet, &stream_tree, graph);

                        // calculate drainage areas
                        drainage_basin.for_each_downstream(|i| {
                            let j = stream_tree.next[i];
                            if j != i {
                                drainage_areas[j] += drainage_areas[i];
                            }
                        });

                        // calculate response times
                        drainage_basin.for_each_upstream(|i| {
                            let j = stream_tree.next[i];
                            let distance: Length = {
                                let (ok, edge) = graph.has_edge(i, j);
                                if ok {
                                    edge
                                } else {
                                    1.0
                                }
                            };
                            let celerity = parameters[i].erodibility * drainage_areas[i].powf(m_exp);
                            response_times[i] += response_times[j] + 1.0 / celerity * distance;
                        });

                        // calculate elevations
                        drainage_basin.for_each_upstream(|i| {
                            let mut new_elevation = elevations[outlet]
                                + parameters[i].uplift_rate
                                    * (response_times[i] - response_times[outlet]).max(0.0);

                            // check if the slope is too steep
                            // if max_slope_func is not set, the slope is not checked
                            if let Some(max_slope) = parameters[i].max_slope {
                                let j = stream_tree.next[i];
                                let distance: Length = {
                                    let (ok, edge) = graph.has_edge(i, j);
                                    if ok {
                                        edge
                                    } else {
                                        1.0
                                    }
                                };
                                let max_slope = max_slope.tan();
                                let slope = (new_elevation - elevations[j]) / distance;
                                if slope > max_slope {
                                    new_elevation = elevations[j] + max_slope * distance;
                                }
                            }

                            changed |= new_elevation != elevations[i];
                            elevations[i] = new_elevation;
                        });
                    });
                }
                else {
                    let above_slopes = (0..num).map(|ia| {
                        let slopes = graph.neighbors_of(ia).iter().filter_map(|ja| {
                            let ediff = elevations[ja.0] - elevations[ia];
                            if ediff > 0.0 {
                                Some((ja.0, (ediff / ja.1).powi(4)))
                            } else {
                                None
                            }
                        }).collect::<Vec<_>>();
                        let slope_sum = slopes.iter().fold(0., |acc, slope| {
                            acc+slope.1
                        });
                        (slopes, slope_sum)
                    }).collect::<Vec<_>>();
                    
                    let below_slopes = (0..num).map(|ia| {
                        let slopes = graph.neighbors_of(ia).iter().filter_map(|ja| {
                            let ediff = elevations[ia] - elevations[ja.0];
                            if ediff > 0.0 {
                                Some((ja.0, (ediff / ja.1).powi(4)))
                            } else {
                                None
                            }
                        }).collect::<Vec<_>>();
                        let slope_sum = slopes.iter().fold(0., |acc, slope| {
                            acc+slope.1
                        });
                        (slopes, slope_sum)
                    }).collect::<Vec<_>>();

                    let mut drainage_areas: Vec<f64> = areas.to_vec();

                    // calculating drainage area
                    for _ in 0..5 {
                        (0..num).for_each(|ia| {
                            let above: &Vec<(usize, f64)> = &above_slopes[ia].0;
                            let area_flown = above.iter().map(|(j, slope)| {
                                if below_slopes[*j].1 > 0.0 {
                                    drainage_areas[*j] * slope / below_slopes[*j].1
                                } else {
                                    0.0
                                }
                            }).sum::<f64>();
                            drainage_areas[ia] = areas[ia] + area_flown; 
                        });
                    }

                    let celerities = (0..num).map(|ia| {
                        parameters[ia].erodibility * drainage_areas[ia].powf(m_exp)
                    }).collect::<Vec<_>>();

                    let mut response_times = vec![0.0; num];

                    for _ in 0..20 {
                        (0..num).for_each(|ia| {
                            let below = &below_slopes[ia].0;
                            let slope_sum = below_slopes[ia].1;
                            let response_time = below.iter().map(|(j, slope)| {
                                response_times[*j] * slope / slope_sum
                            }).sum::<f64>();

                            let distance = below.iter().map(|(j, slope)| {
                                let distance = {
                                    let (ok, edge) = graph.has_edge(ia, *j);
                                    if ok {
                                        edge
                                    } else {
                                        1.0
                                    }
                                };
                                distance * slope / slope_sum
                            }).sum::<f64>();
                            response_times[ia] = response_time + 1.0 / celerities[ia] * distance;
                        });
                    }

                    // calculate elevations
                    (0..num).for_each(|ia| {
                        let new_elevation = elevations[ia]
                            + parameters[ia].uplift_rate * response_times[ia].max(0.0);

                        changed |= new_elevation != elevations[ia];
                        elevations[ia] = new_elevation;
                    });
                }

                // if the elevations of all sites are stable, break
                if !changed {
                    break;
                }
                
                step += 1;
                if let Some(max_iteration) = &self.max_iteration {
                    if step >= *max_iteration {
                        break;
                    }
                }
            }

            elevations
        };

        Ok(model.create_terrain_from_result(&elevations))
    }
}
