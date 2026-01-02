# finite_element_method

![Rust](https://img.shields.io/badge/Rust-stable-orange)
![FEM](https://img.shields.io/badge/Finite%20Element%20Method-FEM-blue)
![Numerical Methods](https://img.shields.io/badge/Numerical%20Methods-assembly%20%26%20stiffness-lightgrey)
![Status](https://img.shields.io/badge/status-learning%20%2F%20building-lightgrey)

A small Rust crate that collects **finite‑element method (FEM) building blocks** used as the
foundation for my browser‑based FEA project [`fea_app`](https://github.com/RomanShushakov/fea_app).

This crate is intentionally **educational and practical**. It focuses on clarity and explicit
data flow (elements → assembly → global system), rather than on being a general‑purpose or
feature‑complete FEM framework.

---

## Overview

The goal of this crate is to make the core ideas of the finite‑element method tangible in code:

- how elements are defined and parameterized,
- how local stiffness contributions are assembled,
- how the global system is prepared for downstream solvers.

It is designed to be read, experimented with, and extended.

---

## What’s inside

- **Mesh and data model**
  - nodes, element connectivity, and basic material / section properties
- **Element implementations**
  - a small set of classic structural elements (e.g. truss / beam / plate variants)
- **Assembly helpers**
  - routines to compute and assemble global stiffness contributions
- **Geometry utilities**
  - small helpers (e.g. planar convex hull) used by element logic

> Solvers live in separate repositories:
>
> - Iterative methods (CG / PCG): [`iterative_solvers`](https://github.com/RomanShushakov/iterative_solvers)
> - Direct solver experiments: [`colsol`](https://github.com/RomanShushakov/colsol)
> - GPU compute & visualization pipeline: implemented in [`fea_app`](https://github.com/RomanShushakov/fea_app)

---

## Design notes

- **Readability first**  
  This code is written to understand FEM internals end‑to‑end, not to hide them behind abstractions.

- **Explicit data flow**  
  Types and functions aim to make mathematical structure and data movement visible.

- **Scoped ambition**  
  This is not a production FEM package. It intentionally omits features such as:
  nonlinear analysis, contact, advanced integration schemes, or a large element catalog.

---

## Documentation

Some background notes and design explanations live in the `docs/` directory:

- [Elements axes direction](./docs/Elements_axes_direction.md)
- [Element types](./docs/Element_types.md)

---

## Relationship to `fea_app`

`fea_app` integrates these building blocks into a complete pipeline:

**mesh input → FEM assembly → solver → visualization (WebGL / WebGPU)**

This crate provides the **FEM side** (elements and assembly).  
Companion crates provide **solvers** and **GPU compute** experiments.

---

## Status

Active learning and exploration.  
The API may evolve as the broader project grows.

---

## License

MIT OR Apache‑2.0 (see `Cargo.toml`).
