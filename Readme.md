
<h1>CAD Design, Meshing, Simulation, Analysis as Code</h1>

##CREATE YOUR VIRTUAL ENVIRONMENT

 -> Install conda on your system

       follow steps at https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

 -> Install everything else 

       conda create -n p_1
       conda activate p_1
       conda install pip
       conda install gmsh pygmsh
       conda install -c conda-forge fenics-dolfinx mpich pyvista
       cd [anaconda3_install_folder]\envs\[p_1]\bin
       before next step, make sure p_1 conda environment is activated
       pip install cadquery
       install cadquery through pip into p_1 conda environment
       cadquery and fenics has some dependecy conflict in conda packages
       cadquery's pip package solves this problem

##WORKFLOW

       Design Modell with cadquery
       Mesh the Modell with Structured (Transfinite) Cubic Mesh
       Simulate Linear Elastic Beam based on https://fenicsproject.org/pub/tutorial/html/._ftut1008.html 
       Analize with Paraview
       Save results in the STEP file
