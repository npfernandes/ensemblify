# The `analysis` module

With the `analysis` module, you can create an interactive graphical dashboard displaying structural information calculated from the conformational ensemble of your protein of interest.

To do this, provide:
- your ensemble in trajectory format;
- your trajectory's corresponding topology file;
- the name you want to use for your protein in the graphical dashboard.

Using the `ensemblify` command in a terminal:

    ensemblify analysis -trj trajectory.xtc -top topology.pdb -tid trajectory_name

Inside a Python script or Jupyter Notebook:

    from ensemblify.analysis import analyze_trajectory
    analyze_trajectory('trajectory.xtc','topology.pdb','trajectory_name')

