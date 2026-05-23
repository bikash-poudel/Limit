# LIMIT

LIMIT (Loosely coupled Modular framework for Isotope Transport) is a modular, host-independent Python framework for stable water isotope transport in soils designed for integration into existing hydrological models.

\---

## Overview

Stable water isotopes are increasingly used as tracers to investigate hydrological processes and water fluxes within the soil–vegetation–atmosphere continuum. However, most existing isotope-enabled hydrological models rely on tightly coupled model architectures in which water, heat, and isotope transport are solved simultaneously within monolithic model structures. Such implementations often limit transferability, reuse, and compatibility across hydrological modeling systems.

LIMIT is designed as a transferable extension layer that enables stable water isotope simulations without requiring intrusive restructuring of existing hydrological models. The framework interacts with host models through the exchange of state variables at configurable temporal resolutions while maintaining consistency with host model dynamics.

The framework supports:

* Liquid-phase isotope transport
* Vapor-phase isotope transport
* Isotopic fractionation processes
* Modular coupling with external hydrological models
* Flexible one-dimensional to three-dimensional model configurations
* Finite-volume based numerical implementation

\---

## Repository Structure

```text
LIMIT/
│
├── README.md
├── LICENSE
├── environment.yml
│
├── src/
│   └── limit/
│       Core implementation of the LIMIT framework
│
└── validation\\\_tests/
    Benchmark and validation experiments presented in the manuscript
```

\---

## Installation

Clone the repository:

```bash
git clone <repository-url>
cd LIMIT
```

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate limit
```

\---

## Quick Start

Run a validation experiment:

```bash
python validation\\\_tests/test\\\_case\\\_1.py
```

The validation experiments reproduce the benchmark and theoretical test cases described in the manuscript.

\---

## Reproducing the Manuscript Results

The validation and benchmark experiments used in the manuscript are provided in the `validation\\\_tests/` directory. Each test case includes the required configuration and execution scripts needed to reproduce the corresponding numerical experiments.

The framework was evaluated using:

* Classical theoretical test cases
* Analytical benchmark comparisons
* Spatio-temporal sensitivity analyses

These experiments demonstrate the numerical consistency and computational performance of the framework.

\---

## Framework Characteristics

LIMIT is designed to:

* Operate independently of a specific hydrological model implementation
* Exchange state variables with external host models at configurable temporal resolutions
* Maintain consistency with host model dynamics while preserving computational efficiency
* Provide optional modules for vapor and heat transfer processes required for physically consistent isotope simulations

The framework is intended as a reusable and extensible isotope transport layer for ecohydrological and land-surface modeling applications.

\---

## Dependencies

All required Python dependencies are listed in:

```text
environment.yml
```

The framework is implemented entirely in Python.

\---

## License

This project is distributed under the MIT License.

\---

## Citation

If you use LIMIT in your research, please cite:

```text
\\\[Your manuscript citation here after publication]
```

\---

## Manuscript

Associated manuscript:

"LIMIT: A modular, host-independent framework for stable water isotope transport in soil"

Submitted to:

Environmental Modelling \& Software



