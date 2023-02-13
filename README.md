## Description

- [Elements axes direction](./docs/Elements_axes_direction.md)
- [Element types](./docs/Element_types.md)
- [Stiffness groups](./docs/Stiffness_groups.md)


## Finite element method

A finite element method module.

### Example

```rust
#[macro_use]
extern crate finite_element_method;

use finite_element_method::{FEM, DOFParameter, ElementForceComponent};

const REL_TOL: f32 = 1e-4;
const ABS_TOL: f32 = 1e-12;


let mut model = FEM::create(REL_TOL, ABS_TOL, 10);

model.add_node(1, 0.0, 0.0, 0.0)?;
model.add_node(2, 0.0, 30.0, 0.0)?;

model.add_truss(1, 1, 2, 1e6, 2.0, None)?;

model.add_displacement(1, DOFParameter::X, 0.0)?;

model.add_concentrated_load(2, DOFParameter::Y, 100.0)?;

let separated_stiffness_matrix = model.separate_stiffness_matrix()?;
let r_a_vector = model.compose_r_a_vector(separated_stiffness_matrix.get_k_aa_indexes())?;
let u_b_vector = model.compose_u_b_vector(separated_stiffness_matrix.get_k_bb_indexes())?;

let u_a_vector = model.find_ua_vector(
    &separated_stiffness_matrix, &r_a_vector, &u_b_vector,
)?;
let r_r_vector = model.find_r_r_vector(
    &separated_stiffness_matrix, &u_a_vector, &u_b_vector,
)?;

model.compose_global_analysis_result(
    separated_stiffness_matrix.get_k_aa_indexes(), 
    separated_stiffness_matrix.get_k_bb_indexes(), 
    &u_a_vector,
    &r_r_vector,
)?;

let mut global_analysis_result = model.extract_global_analysis_result()?;
global_analysis_result.sort_by(
    |(n_1, dof_1, _, _), (n_2, dof_2, _, _)| 
    (n_1, dof_1).partial_cmp(&(n_2, dof_2)).unwrap()
);
let global_analysis_result_expected = vec![
    (1, DOFParameter::X, 0.0, -100.0),
    (1, DOFParameter::Y, 0.0, 0.0),
    (1, DOFParameter::Z, 0.0, 0.0),
    (1, DOFParameter::ThX, 0.0, 0.0),
    (1, DOFParameter::ThY, 0.0, 0.0),
    (1, DOFParameter::ThZ, 0.0, 0.0),
    (2, DOFParameter::X, 0.0014999999, 100.0),
    (2, DOFParameter::Y, 0.0, 0.0),
    (2, DOFParameter::Z, 0.0, 0.0),
    (2, DOFParameter::ThX, 0.0, 0.0),
    (2, DOFParameter::ThY, 0.0, 0.0),
    (2, DOFParameter::ThZ, 0.0, 0.0),
];

let elements_analysis_result = model.extract_elements_analysis_result()?;
let elements_analysis_result_expected = vec![
    (1, vec![(ElementForceComponent::ForceR, 100.0)]),
];

assert_eq!(global_analysis_result, global_analysis_result_expected);
assert_eq!(elements_analysis_result, elements_analysis_result_expected);
```
