
<h1>CAD Design, Meshing, Simulation, Analysis as Code</h1>

<h2>CREATE YOUR VIRTUAL ENVIRONMENT</h2>

Install conda on your system

- Follow steps at: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Install everything else 

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
       cd [to the project folder]

<h2>WORKFLOW</h2>

1. Design Modell with cadquery
2. Mesh the Modell with Structured (Transfinite) Cubic Mesh
3. Simulate Linear Elastic Beam based on https://fenicsproject.org/pub/tutorial/html/._ftut1008.html 
4. Analize with Paraview
5. Save results in the STEP file
