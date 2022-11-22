use std::{cmp::min, collections::VecDeque};

// row major
// 0,0 bottom left
// things are only supported on the ground
const WIDTH: usize = 5;
const HEIGHT: usize = 5;
const WEIGHT: f32 = 10.0;
type Grid<A> = [[Option<A>; WIDTH]; HEIGHT];

fn display_grid<T>(grid: &Grid<T>, print: impl Fn(&T) -> char) -> String {
    let mut out = String::new();
    out.push('+');
    for _ in 0..WIDTH {
        out.push('-');
    }
    out.push_str("+\n");

    for row in grid.iter().rev() {
        out.push('|');
        for cell in row {
            out.push(if let Some(content) = cell {
                print(content)
            } else {
                ' '
            });
        }
        out.push_str("|\n");
    }
    out.push('+');
    for _ in 0..WIDTH {
        out.push('-');
    }
    out.push_str("+\n");

    out
}

fn display_unit_grid(grid: &Grid<()>) -> String {
    display_grid(grid, |_| 'X')
}
fn display_int_grid(grid: &Grid<u32>) -> String {
    display_grid(grid, |content| {
        // simply don't have a depth greater than 9 when using this visualisation
        format!("{}", content).chars().next().unwrap()
    })
}

fn convert_to_distance_from_ground(grid: &Grid<()>) -> Grid<u32> {
    // (x,y)
    let mut queue = VecDeque::new();
    // if None, then not visited yet. If Some, then best estimate seen so far.
    let mut out = [[None; WIDTH]; HEIGHT];

    // add all ground nodes to the front of the queue
    for x in 0..WIDTH {
        queue.push_back((x, 0));
    }

    // it's your friend: BFS
    while let Some(ix) = queue.pop_front() {
        // no block = nothing to do
        if index(grid, ix).is_some() {
            // work out our best distance to ground
            // println!("({}, {})", x,y);
            let dist_from_ground =
                // if on ground, it's 0
                if bottom_edge(ix) {
                    0
                } else {
                    // candidates in 4 directions
                    neighbour_indices(ix)
                        .into_iter()
                        .filter_map(|ix| index(&out, ix).as_ref())
                        .min()
                        .expect("no neighbours??") + 1
                };
            *index_mut(&mut out, ix) = Some(dist_from_ground);

            // println!("({}, {}): {}", x, y, dist_from_ground);

            // if a neighbour has been visted before, update its shortest path if needed
            // if not, set its shortest path and add it to the queue
            neighbour_indices(ix)
                .into_iter()
                .filter(|&ix| index(grid, ix).is_some())
                .for_each(|ix| match *index(&out, ix) {
                    Some(neighbour_distance) => {
                        *index_mut(&mut out, ix) =
                            Some(min(neighbour_distance, dist_from_ground + 1))
                    }
                    None => {
                        *index_mut(&mut out, ix) = Some(dist_from_ground + 1);
                        queue.push_back(ix);
                    }
                });
        }
    }

    out
}

fn order_by_rev_distance_from_ground(grid: &Grid<u32>) -> Vec<(usize, usize)> {
    let mut out = Vec::new();

    for (y, row) in grid.iter().enumerate() {
        for (x, cell) in row.iter().enumerate() {
            if cell.is_some() {
                out.push((x, y))
            }
        }
    }

    // FIXME: we end up looking up this info a bunch anyway why throw it away here
    out.sort_unstable_by_key(|&ix| index(grid, ix).unwrap());
    out.reverse();
    out
}

// FIXME: gotta care about center of mass of the neighbours when calculating torque, rip
#[derive(Copy, Clone, Debug, Default)]
struct Force {
    linear: (f32, f32), // 2d vector
    torque: f32,        // pos is clockwise, about center
    supporting_neighbours: u32,
}

// input = grid, list of coords sorted by reverse distance from ground
fn calculate_moments(grid: &Grid<u32>, blocks: &Vec<(usize, usize)>) -> Grid<Force> {
    let mut out: Grid<Force> = [[None; WIDTH]; HEIGHT];
    for &ix in blocks {
        // println!("calculating ({},{})", x, y);
        let mut self_force = index(&out, ix).unwrap_or_default();
        self_force.linear.1 -= WEIGHT;
        let self_dist = index(grid, ix).unwrap();

        // okay so we look at neighbours.
        // any further than us from the ground we assume are already calculated
        let neighbours = neighbour_indices(ix)
            .into_iter()
            .filter_map(|ix| index(grid, ix).map(|dist| (ix, dist)));

        let mut dependent_neighbour_count = 0;
        for (neighbour_ix, dist) in neighbours.clone() {
            // println!("looking at neighbour ({}, {}): {:?}", x, y, out[y][x]);
            if dist <= self_dist {
                dependent_neighbour_count += 1;
            } else {
                let neighbour_force = index(&out, neighbour_ix).unwrap();
                // fully calculated, work out the torque it applies to us
                if ix.0 == neighbour_ix.0 {
                    // it's above or below, so only care about x dir torque
                    if neighbour_ix.1 < ix.1 {
                        // it's below, so right is anticlockwise
                        // FIXME: should this be half supporting neighbours because each one can only hold on with the half
                        //        of it that's in tension
                        self_force.torque -= neighbour_force.linear.0
                            / (neighbour_force.supporting_neighbours as f32);
                    } else {
                        // it's above, so right is clockwise
                        self_force.torque += neighbour_force.linear.0
                            / (neighbour_force.supporting_neighbours as f32);
                    }
                } else {
                    // it's left or right, so only care about y dir torque
                    if neighbour_ix.0 < ix.0 {
                        // it's left, so upwards is clockwise
                        self_force.torque += neighbour_force.linear.1
                            / (neighbour_force.supporting_neighbours as f32);
                    } else {
                        // it's right, so upwards is anticlockwise
                        self_force.torque -= neighbour_force.linear.1
                            / (neighbour_force.supporting_neighbours as f32);
                    }
                }
            }
        }
        self_force.supporting_neighbours = dependent_neighbour_count;
        *index_mut(&mut out, ix) = Some(self_force);

        for (ix, dist) in neighbours {
            if dist <= self_dist {
                // not yet calculated, add our forces
                if index(&out, ix).is_none() {
                    *index_mut(&mut out, ix) = Some(Force::default());
                }
                let neighbour_force = index_mut(&mut out, ix).as_mut().unwrap();
                (*neighbour_force).linear.0 +=
                    self_force.linear.0 / (dependent_neighbour_count as f32);
                (*neighbour_force).linear.1 +=
                    self_force.linear.1 / (dependent_neighbour_count as f32);
            }
        }
    }

    out
}

fn left_edge((x, _): (usize, usize)) -> bool {
    x == 0
}
fn right_edge((x, _): (usize, usize)) -> bool {
    x >= WIDTH - 1
}
fn bottom_edge((_, y): (usize, usize)) -> bool {
    y == 0
}
fn top_edge((_, y): (usize, usize)) -> bool {
    y >= HEIGHT - 1
}
fn neighbour_indices((x, y): (usize, usize)) -> Vec<(usize, usize)> {
    let mut out = vec![];
    if !left_edge((x, y)) {
        out.push((x - 1, y))
    }
    if !right_edge((x, y)) {
        out.push((x + 1, y))
    }
    if !bottom_edge((x, y)) {
        out.push((x, y - 1))
    }
    if !top_edge((x, y)) {
        out.push((x, y + 1))
    }
    out
}

// just because I know I will fuck this up one day. maybe I should have done column major
fn index<T>(grid: &Grid<T>, (x, y): (usize, usize)) -> &Option<T> {
    &grid[y][x]
}
fn index_mut<T>(grid: &mut Grid<T>, (x, y): (usize, usize)) -> &mut Option<T> {
    &mut grid[y][x]
}
fn demo(test_grid: Grid<()>) {
    println!("{}", display_unit_grid(&test_grid));
    let int_grid = convert_to_distance_from_ground(&test_grid);
    // println!("{:?}", int_grid);
    println!("{}", display_int_grid(&int_grid));

    let sorted_grid = order_by_rev_distance_from_ground(&int_grid);
    println!("{:?}", sorted_grid);

    let forces = calculate_moments(&int_grid, &sorted_grid);
    println!("{:?}", forces);
}
fn main() {
    // this is upside down lol
    let test_grid = [
        [None, None, Some(()), None, None],
        [Some(()), Some(()), Some(()), None, None],
        [None, None, None, None, None],
        [None, None, None, None, None],
        [None, None, None, None, None],
    ];
    demo(test_grid);

    let test_grid = [
        [Some(()), None, Some(()), None, None],
        [Some(()), Some(()), Some(()), None, None],
        [None, Some(()), None, Some(()), None],
        [None, Some(()), Some(()), Some(()), None],
        [None, None, None, None, None],
    ];
    demo(test_grid);

    let test_grid = [
        [Some(()), None, None, None, Some(())],
        [Some(()), None, None, None, Some(())],
        [Some(()), Some(()), Some(()), Some(()), Some(())],
        [None, None, None, None, None],
        [None, None, None, None, None],
    ];
    demo(test_grid);

    let test_grid = [
        [Some(()), None, None, None, Some(())],
        [Some(()), Some(()), None, Some(()), Some(())],
        [Some(()), Some(()), Some(()), Some(()), Some(())],
        [None, None, None, None, None],
        [None, None, None, None, None],
    ];
    demo(test_grid);
}
