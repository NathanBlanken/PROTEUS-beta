# Contrast-enhanced ultrasound simulator for microbubbles in vascular flow
- Collaboration between University of Twente (Physics of Fluids, Engineering Fluid Dynamics) and Delft University of Technology (Imaging Physics)

## Installation of the simulator

Download (`code` > `Download ZIP` > unzip) or clone the repository to your device.

Download the k-Wave binary files from (www.k-wave.org/download.php). Project developers can also download the binaries from SURFdrive. Copy the binary files into the folder `k-wave-toolbox-version-1.3/k-Wave/binaries`. For custom installation paths, see _Custom installation_ below.

The flow simulation is saved in VTU format, which is converted to a MATLAB file named `vtu.mat`. For the renal tree, download the file `vtu.mat` from SURFdrive and copy it to the folder `geometry_data/renal_tree`. More geometries will be added in the future.


## Flow solver module

This is the general link: https://github.com/apes-suite.

The link to LBM solver documentation specifically is: https://geb.inf.tu-dresden.de/doxy/musubi/index.html.

## Running the simulator on a local device

To open the GUI type `MainGUI` in the MATLAB command window.

Modify the simulation settings by opening the relevant pop up windows (Microbubbles, Simulation Parameters, Geometry and streamlines, Medium, Transducer, Transmit, Acquisition). 

Run the simulation, either pressure field simulation or RF data simulation. The results (pressure maps or RF data) will be stored in `RESULTS/<timestamp>`. The ground truth data used for the simulation is stored in `ground_truth_frames/<timestamp>`. Upon first run, the folder `RESULTS` and `ground_truth_frames` are automatically created.


## Running the acoustic module remotely

To run a simulation on a remote system without graphical interface, copy the repository to the remote system and perform the steps below.

Pressure field simulation:
- Run `MainGUI` on your local device to set the simulation settings. Save the settings by clicking `Save settings`.
- Copy the saved settings file to the remote system.
- On the remote system, add the acoustic module to the MATLAB path with the command `addpath('acoustic-module')`.
- Run the function `main_pressure_field`. See the function's help (`help main_pressure_field`) for help on the input arguments. 

RF data simulation:
- Copy the repository to the remote system.
- Run `MainGUI` on your local device to set the simulation settings. Save the settings by clicking `Save settings`.
- Precompute ground truth data (see _Precomputing ground truth data_ below).
- Copy the ground truth data and the saved settings file to the remote system.
- On the remote system, add the acoustic module to the MATLAB path with the command `addpath('acoustic-module')`.
- Run the function `main_RF`. See the function's help (`help main_RF`) for help on the input arguments. 


## Precomputing ground truth data
To precompute ground truth microbubble positions (i.e. prior to acoustic simulation), perform the following steps):
- Run `MainGUI`.
- Open the Acquisition window by clicking `Modify`.
- Tick the box `Reuse or precompute ground truth positions`.
- Click `Precompute ground truth`.

## Custom installation (not recommended)
When installing components (such as the k-Wave binaries) in different locations, you need to modify the file `PATHS.mat`:
- Open the file `PAHTS.mat`.
- Enter `'full path to module'` in the relevant `Path` field. Do NOT modify the `Description` field.
- Save the table as `PATHS.mat`.

