# main.py

import cadquery
import gmsh
import dolfinx
import sys

from cqplugin.export_step import export_step_ap214

def main() -> int:

    # Scaled variables
    L = 1.0; W = 0.2
    mu = 1
    rho = 1
    delta = W/L
    gamma = 0.4*delta**2
    beta = 1.25
    lambda_ = beta
    g = gamma

    # create cadquery model
    cad_model = cadquery.Workplane("xy").box(L, W, W).toOCC()

    # mesh it 
    

    # simulation
    # analysis
    # saving

    return 0 

if __name__ == '__main__':
    sys.exit(main())