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
    Number of threads. Default=1. Threads do not help for Riemann since it is algorithmic and deterministic.

-b [NUM]
    Batch size. Default=1e7. Each thread will run 'b' batches.

-a [NUM]
    Start of integral. Default = -10.

-z [NUM]
    Endpoint of integral. Default = 10.

-f [FUNCTION NAME]
    {identity | normal | normaln} Name of test function to integrate. Default = identity.

-i [INTEGRATOR]
    { riem | stoc } Name of integration algorithm to use. Default = riemann.


Example, integrating Normal from 0 to 2 with Rieman, batch = 1e6:

`./bin/integrator -b 1e6 -f normal -a 0 -z 2 -i riem`

Example, integrating Normal from 0 to 2 with Stochastic, batch = 1e7, 4 threads:

`./bin/integrator -b 1e7 -f normal -a 0 -z 2 -i stoc -t 4`

Some results comparing riemann and stochastic:
```
Normal, -2 to 2 (Ans = 0.95449973610364158, Python Scipy normal.cdf)
Batch = 1e8
Riemann:    0.954499736104 (perfect to 12 dp)
Stochastic: 0.954539838414 (difference 4.20e-03%)

Batch = 1e3
Riemann:    0.95449988007 (difference 1.508e-5%)
Stochastic: 0.94922105965 (difference 0.553%)

Batch = 1e2
Riemann:    0.9545141330226 (difference -0.0015%)
Stochastic: 1.0090704489544 (difference -5.717%)
```

Conclusion: Riemann integration is much faster to execute (even with threading) and also much more accurate for a
given batch size as compared to stochastic integration (naive implementation).

