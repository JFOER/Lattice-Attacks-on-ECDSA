# Lattice-Attacks-on-ECDSA
This repository encompasses the practical realization of Lattice Attacks against the Elliptic Curve Digital Signature Algorithm (ECDSA), as outlined in the following research paper:
"Back to Construct a More Effective Lattice: a Generic Attack Model on ECDSA"

# How to Run the script
1. Directly copy the code into Pycharm(or other IDE). Please note that our code requires the support of the fpylll package. Ensure that it is properly installed for optimal functionality.
2. Or directly copy the code into the sage terminal.

# data.7z
This compressed file contains 100 folders, with each folder housing the signature data collected from 1000 signatures using the same private key, including values such as 'r', 's', etc. The private keys differ among different folders. Please ensure that the file path is correctly set up before use.

# Attack MSBs or LSBs（variant 1）.py
This script implements Variant 1 to attack the MSB or LSB leakage models, where the number of signatures, the number of leaked bits, and other parameters can be specified.
