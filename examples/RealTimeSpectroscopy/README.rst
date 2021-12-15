Examples - RealTimeSpectroscopy
===============================

The binary: RTSpectroscopyDriver
--------------------------------

The exploration along trajectories of concatenated XYZ-coordinates is
possible with the RTSpectroscopyDriver binary. This is a thin wrapper
around Sparrow parsing the structures in the trajectories and 
carrying out calculations on them sequentially.

Two ingredients are necessary to carry out the exploration along a trajectory, 
the trajectory and the input file.

The trajectory: example_trajectory.xyz
--------------------------------------

A trajectory is defined by a concatenated sequence of xyz coordinate files. 
As an example of the format, we provide ``example_trajectory.xyz``.

The input file: input.yaml
--------------------------

The program must have some information about how the calculation must be
carried out. In particular, the desired properties need to be specified.
This is done in the ``input.yaml`` file. 

The exploration
---------------

The exploration is then started by calling

```
RTSpectroscopyDriver input.yaml
```
