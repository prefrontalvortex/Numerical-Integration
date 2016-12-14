# INSTALLING
## Prerequisites:
You must have CMake v 3.6 or better and GCC 6 or higher. I can't guarantee it'll work otherwise. Linux is deal but
this SHOULD work on Mac.

Unzip the archive and cd into the project directory.
```
mkdir build
cmake ..
make
```

cd back to the main directory and run with

`./bin/integrator [OPTIONS]`
`./bin/integrator [-t N_THREADS] [-b BATCHSIZE] [-a START] [-z END] [-f {identity | normal | normaln}]
                    [-i { riem | stoc }]`


-t [INT]
    Number of threads. Threads do not help for Riemann since it is algorithmic and deterministic.

-b [NUM]
    Batch size. Each thread will run 'b' batches.

-a [NUM]
    Start of integral

-z [NUM]
    Endpoint of integral

-f [FUNCTION NAME]
    {identity | normal | normaln} Name of test function to integrate

-i [INTEGRATOR]
    { riem | stoc } Name of integration algorithm to use.
