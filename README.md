# QuantSatSim_v1.0.0
Modelling for different SatQKD systems based on the SatQuMA v2.0.0 software developed and released by the Space Quantum Technologies team in the Computational Nonlinear & Quantum Optics group in the Department of Physics, University of Strathclyde.

This software is designed to simulate a SatQKD communication using the simplified BB84 and was developed by the Quantum Photonics Laboratory (QuLab) in the Instituto de Telecomunicações in the Department of Electrical Engineering and Computers, Instituto Superior Técnico, University of Lisbon. The software is part of a project described in the article [1]

[1] Mendes P.N., Teixeira G.L., et al. "Optical payload design for downlink quantum key distribution and keyless communication using CubeSats." EPJ Quantum Technology 11.1 (2024): 1-19.

# Getting Started
## Setup
The necessary Python packages are stated in a `requirements.txt` file. To install them just run 
```bash
pip install -r requirements.txt
```

## Run the Code
Edit the input parameters in `inputfiles/input.txt` before running the code.
To run the simulations start with the file `QuantSatSim.py`.
