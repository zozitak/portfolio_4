# main.py

import cadquery as cad
import gmsh
import dolfinx
import sys

from cqplugin.export_step import export_step_ap214

def main() -> int:

    # Scaled variables
    L = 1; W = 0.2
    mu = 1
    rho = 1
    delta = W/L
    gamma = 0.4*delta**2
    beta = 1.25
    lambda_ = beta
    g = gamma

    # create cq model

    # mesh it 
    # simulation
    # analysis
    # saving
    return 0 

if __name__ == '__main__':
    sys.exit(main())