# num_methods

[![Matlab Version](https://img.shields.io/badge/MATLAB-R2024%2B-blue)](https://www.mathworks.com/products/matlab.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md)

A repository of academic Matlab scripts ranging various areas of numerical methods. It aims to provide clean and understandable code, more tailored towards academic applications rather than production code.

## Scope and Audience

This repository is mainly intended for undergraduate and early graduate students in numerical methods, scientific computing, and engineering courses, as well as, anyone looking to understand how numerical algorithms work internally. The focus is on clarity and educational value rather than performance or production readiness.

It covers the basic methods any engineer or mathematician should know (Newton, RK4, LU factorization,...), as well as some more exotic methods (AX19, HMT,...). It can also aid in understanding the patterns and procedures to extend methods into iterative schemes, as many of the basic programs are also adapted into iterative variations.

Finally, the example files document the programs used, showing examples for every program, as well as some of the limitations inherent to these methods.

## Design Philosophy

This repository aims to have

- Readable implementations over optimized code
- Minimal dependencies
- Consistent function interfaces where possible
- Input validation on all function
- Well-documented functions using MATLAB docstrings

The goal is to make each method easy to follow and adapt for coursework or study.

## Contents

Currently, this repository has folders dedicated to

- Root Finding
- Integration
- Interpolation
- Linear Systems
- Non Linear Systems
- ODE

See the [INDEX](Docs/INDEX.md) file for the whole list of methods, as well as their documentation.

## Quick start / Setup

In order to use this programs, you have to

1. Clone or download the repository
2. Add it to your MATLAB path
3. Verify installation (Run `startup` and `help Newton`, for example)

The folder contains a `startup.m` file, so running `startup` in the terminal will add all subfolders to your working directory until the program is closed. This way, programs can be run without their folder having to be open.

## Usage notes

Most programs work on base MatLab, but certain ones (gaussian cuadrature, interpolation) require the symbolic math toolbox. All programs are documented with docstrings, so typing `help foo` for any one of the program's names will show their syntax and arguments. Make sure to have the files in the current MatLab working directory (or one added to the path) for them to work.

## Limitations

The code in this repository

- Is not optimized for large-scale or high-performance computations
- Assumes mathematically well-posed problems (disregards effects of numerical error)

## License

This project is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).

See the [LICENSE](LICENSE.md) file for details.

If you use this repository, please cite following the [citation](citation.cff).
