   This program calculates the Mean Squared Displacement (MSD) of selected
   atoms (e.g., Li) from a LAMMPS trajectory file. It reads configuration
   parameters from an input file, processes the trajectory, and outputs the
   MSD as a function of time.

 Main Features:
   - Reads simulation and analysis parameters from 'inp_msd_tot.dat'.
   - Reads atom coordinates from a LAMMPS trajectory file.
   - Selects atoms of a specified type (e.g., lipids) for MSD calculation.
   - Computes MSD in x, y, z, xy, and total directions over time.
   - Outputs results to a specified file. Input Files:
   - inp_msd_tot.dat: Configuration file containing paths, filenames, and
     analysis parameters (e.g., atom types, distance intervals).
   - LAMMPS trajectory file: Contains atomic coordinates for each frame.
   -  Output Files:
   - MSD results written to the file specified by 'output_file' in the config.

 Key Variables:
   - na: 3D array storing coordinates of selected atoms over all frames.
   - msd_x, msd_y, msd_z, msd_xy, msd_tot: Arrays storing MSD values.
   - li_id: Maps atom IDs to indices for selected atom type.
   - Parameters such as 'lip', 'lig', 'dist_int', 'r_iz', 'l_iz' are read
     from the configuration file.

Usage:
   1. Prepare 'inp_msd_tot.dat' with correct paths and parameters.
   2. Place the trajectory file in the specified location.
   3. Compile and run the program.
   4. MSD results will be written to the output file.

 Notes:
   - The program assumes a specific LAMMPS trajectory format.
   - The number of frames and atoms must be consistent with the input files.
   - The program is tailored for a fixed number of frames and atoms, as set
     by the parameters at the top of the code.
