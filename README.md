
<h1>CAD Design, Meshing, Simulation, Analysis as Code</h1>

<h2>Introduction</h2>

Demonstration of a text-based engineering workflow aimed at showcasing proper version control practices in engineering.

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

<h3>1. Design Modell with cadquery</h3>

<img src="/assets/cad.png">

<h3>2. Mesh the Modell with Structured (Transfinite) Cubic Mesh</h3>

<img src="/assets/gmsh.png">

<h3>3. Simulate Linear Elastic Beam based on https://fenicsproject.org/pub/tutorial/html/._ftut1008.html </h3>

       Scenario:

       A cantilever beam is clamped at one end

                  .+------------------------+
                .' |                      .'|
               +---+--------------------+'  |  ↓ gravity
       clamped |   |                    |   |
               |  ,+--------------------+---+
               |.'                      | .'
               +------------------------+'

       It is subject to the load due to its own weight and will
       deflect accordingly. Under an assumpation of small
       deformation the material follows linear elasticity.

<h3>4. Analize with Paraview</h3>

<img src="/assets/paraview.png">
Maximum Displacement: 0.23813392114555948 unit
</br>
<h3>5. Save results in the STEP AP214 standard file format.</h3>

Test: test_step_ap214_export ... FAIL

Sadly based on Opencascade documentation (https://dev.opencascade.org/doc/occt-7.6.0/overview/html/occt_user_guides__step.html#autotoc_md254)</br>
Opencascade doesnt export notes or custom data in any of the STEP formats. 

For this purpose in the future either OCCT should be echanced or another CAD Kernel could be used.</br>
This Feature would be a great asset of future CAD Kernels. 