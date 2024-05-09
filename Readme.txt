NAME

CAD Design, Meshing, Simulation, Analysis as Code

DESCRIPTION

State of the art CAD design softwares don't come with version control system. 
This aproach attempts to solve this problem on a maintainable way. 
Additionaly gives the chance to use the STEP AP214 standard file format, 
to pack extra manufacturing & simulation data into CAD models.

CREATE YOUR VIRTUAL ENVIRONMENT

 -> Install conda on your system
        follow steps at https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

 -> Install everything else 

        conda create -n p_1
        conda activate p_1
        conda install pip
        conda install gmsh pygmsh
        conda install -c conda-forge fenics-dolfinx mpich pyvista
        cd [anaconda3_install_folder]\envs\[p_1]\bin
        #before next step, make sure p_1 conda environment is activated
        pip install cadquery
        #install cadquery through pip into p_1 conda environment
        #cadquery and fenics has some dependecy conflict in conda packages
        #cadquery's pip package solves this problem

