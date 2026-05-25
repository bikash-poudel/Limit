# LIMIT

LIMIT (Loosely coupled Modular framework for Isotope Transport) is a modular, host-independent Python framework for stable water isotope transport in soils designed for integration into existing hydrological modeling systems.

The framework enables isotope transport simulations without requiring tightly coupled monolithic implementations in which water, heat, and isotope transport are solved within a single model architecture.

---

## Overview

Stable water isotopes are widely used as tracers for investigating hydrological processes and water fluxes within the soil–vegetation–atmosphere continuum. Despite their scientific importance, isotope-enabled hydrological models remain relatively limited because most existing implementations rely on tightly coupled model architectures that are difficult to transfer, reuse, and extend.

LIMIT addresses this limitation through a loosely coupled and modular design philosophy. The framework is intended to function as a reusable isotope transport extension layer that interacts with external hydrological host models through the exchange of state variables and fluxes at configurable temporal resolutions.

The framework supports:

- Liquid-phase isotope transport
- Vapor-phase isotope transport
- Isotopic fractionation processes
- Finite-volume based transport simulations
- Flexible one-dimensional to three-dimensional model domains
- Modular integration with external hydrological models

---

## Host Model Dependency

LIMIT is designed to operate in conjunction with an external hydrological host model. The framework relies on the host model to provide hydrological states and fluxes required for isotope transport calculations.

The current validation examples use the Catchment Modelling Framework (CMF) as an example host-model implementation. However, LIMIT is designed to remain independent of a specific hydrological model architecture and can be adapted to other hydrological modeling systems.

CMF is not included as a core dependency of the framework and must be installed separately to reproduce the CMF-based validation examples.

---

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
└── validation_tests/
    CMF-based benchmark and validation experiments
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/bikash-poudel/Limit.git
cd LIMIT
```

Create the conda environment:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate limit
```

Install LIMIT in editable mode:

```bash
pip install -e .
```

If using Git Bash on Windows:

```bash
source activate limit
```

---

## Environment Management

The conda environment specification is provided in:

```text
environment.yml
```

This file contains the core Python dependencies required for the LIMIT framework.



## Running Validation Examples

All commands assume execution from the repository root directory.

Example:

```bash
python validation_tests/iso_cmf.py
```

The validation examples reproduce the benchmark and theoretical test cases described in the associated manuscript.

Some validation examples require external host-model dependencies such as CMF. These dependencies are not included in the core LIMIT environment and must be installed separately.

---

## Reproducing the Manuscript Results

The benchmark and validation experiments presented in the manuscript are located in:

```text
validation_tests/
```

These experiments include:

- Classical theoretical test cases
- Analytical benchmark tests

The provided examples are intended to reproduce the numerical behavior and validation results discussed in the manuscript.

---

## Code Structure and Implementation Guidelines

The framework follows a modular implementation philosophy intended to simplify integration with different hydrological host models.

General implementation guidelines:

- Keep isotope transport components independent from host-model-specific code
- Separate numerical transport logic from host-model coupling interfaces
- Avoid hard-coded dependencies on a particular hydrological model structure
- Preserve modularity and extensibility of transport components
- Maintain compatibility with finite-volume based transport formulations

When extending the framework:

- Clearly separate framework functionality from validation scripts
- Document newly introduced dependencies
- Preserve reproducibility of benchmark examples
- Maintain consistency with the framework's host-independent design philosophy

---

## Dependencies

The framework is implemented entirely in Python.

Core dependencies are listed in:

```text
environment.yml
```

Additional dependencies may be required for specific host-model validation examples.

---

## License

This project is distributed under the MIT License.

---

## Citation

If you use LIMIT in your research, please cite:

```text

A DOI and formal citation information will be added after publication.

---

## Associated Manuscript

"LIMIT: A modular, host-independent framework for stable water isotope transport in soil"

Submitted to:

Environmental Modelling & Software
