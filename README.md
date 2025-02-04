# DD2444-Project: A Distributed Computing Implementation of a 3-dimensional Extension of the Vicsek Model
Project (graded, 7.5 credits) for the Project Course in Scientific Computing (DD2444).
More information about all the steps can be found in the project report.

## Images

A flock of many aligned particles.
![A flock of many aligned particles](https://i.imgur.com/JtO8Dpb.png)

A sparser area of particles with higher concentrations surrounding. Particles that are further away from each other have a higher chance of going in a different direction. 
![A sparser area of particles with higher concentrations surrounding](https://i.imgur.com/0sWmwmU.png)

## Files
- **proj_common.h**
    - Contains all functions and constants (parameters) neccessary to all of the different implementations of the Vicsek model.
- **proj_base.c**
    - Contains a simple implementation of the project, fully functioning but without optimizations.
- **proj_mpi_io.c**
    - Contains an implementation utilizing MPI IO to improve the speed of writing data to file.
- **proj_boxes.c**
    - Contains an implementation using both MPI IO as well as cutting up the entire volume into "boxes" of at least the size of the interaction radius.
- **Godot Visualization/**
    - Folder that contains project files for the Godot visualization of the results.  

## Running Locally

This assumes you are running **Ubuntu** for your local computer. 

Clone the project

```bash
  git clone git@github.com:Teonord/DD2444-Project.git
```

Go to the project directory

```bash
  cd DD2444-Project
```

Change parameters like number of birds and timestep in **proj_common.h**

```bash
  emacs proj_common.h
```

### Running Code

Install MPICH

```bash
  sudo apt install mpich
```

Compile the files

```bash
  mpicc -o <filename>.out <filename>.c -fopenmp -lm
```

Set amount of OpenMP threads (if needed)

```bash
  export OMP_NUM_THREADS=<threads>
```

Run code with number of processes

```bash
  mpiexec -n <processes> ./<filename>.out > res
```
## Running on Dardel
Follow guide in Running Locally to download and set up code.

Copy code to Dardel

```bash
  scp <DD2444-Project> <username>@dardel.pdc.kth.se:Private/.
```

After entering Dardel and moving to proper folder, compile code.

```bash
   cc -o <filename>.out <filename>.c -fopenmp -lm
```

Create Job Script for running code.

```bash
   touch <filename>.sh
   emacs <filename>.sh
```

Enter information neccessary for your batched run of the code.

```bash
    #!/bin/bash -l

    #SBATCH -A <allocation_code>

    #SBATCH -J <filename>

    #SBATCH -p main

    #SBATCH -n <n>

    #SBATCH --nodes=<n>

    #SBATCH --ntasks-per-node=<n>

    #SBATCH --cpus-per-task=<n>

    #SBATCH -t 00:05:00

    export OMP_NUM_THREADS=<threads>
    export OMP_PLACES=cores

    srun ./<filename>.out > ./Results/res_<filename>
```

Save and exit using Ctrl-x, Ctrl-c, y, Enter.

Run the batched Job and wait for run to be completed.

```bash
   sbatch <filename>.sh
```

To check progress, you can check the queue.

```bash
   squeue -u <username>
```


## Visualization
Load the Godot Visualization folder as a project in Godot 4. 

Put the Results file and the Output.dat file in the project folder. 

Select the Controller node and change properties as desired. If using the **proj_base.c** code, enable "Skip Other" and disable "Pure Doubles". If using any other implementation, enable "Pure Doubles" and disable "Skip Other".

Run the Main Scene. A visualization controllable by mouse, WASD, shift, and space will now appear. 
