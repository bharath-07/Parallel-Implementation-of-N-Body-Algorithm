# Project Readme - N-Body Problem
This README contains the instructions to run all the algorithms that we have
developed for the N-body problem.

We tested our code on the crunchy machines in CIMS, NYU. More specifically, 
we used crunchy5 for our project.


## Pre-requisites
1. Before running any of the codes, load gcc-9.2 : ```module load gcc-9.2```.


## Generating inputs
1. Compile the input generation code as follows : ```g++ gen_nbody.cpp```

2. To generate 1000 bodies, please run : ```./a.out 1000```

3. This generates 1000 randomly placed bodies in a 3D space with random masses in a file 1000.txt.


## Sequential Version
1. Compile the input generation code as follows : ```g++ serial-nbody2.cpp```

2. To run the sequential code on 1000 bodies, please run : ```./a.out 1000.txt```


## Parallel Version - 1
1. Compile the input generation code as follows : ```g++ -fopenmp parallel_1.cpp```

2. To generate 1000 bodies, using 32 parallel threads please run : ```./a.out 1000.txt 32```


## Parallel Version - 2
1. Compile the input generation code as follows : ```g++ -fopenmp parallel_2.cpp```

2. To generate 1000 bodies, using 32 parallel threads please run : ```./a.out 1000.txt 32```


## Parallel Version - 3
1. Compile the input generation code as follows : ```g++ -fopenmp parallel_3.cpp```

2. To generate 1000 bodies, using 32 parallel threads, and if you want to
   divide the n bodies into 8 chunks, please run : ```./a.out 1000.txt 32 8```


## Parallel Version - 4
1. Compile the input generation code as follows : ```g++ -fopenmp parallel_4.cpp```

2. To generate 1000 bodies, using 32 parallel threads, and if you want to
   divide the n bodies into 8 chunks, please run : ```./a.out 1000.txt 32 8```


## Parallel Version - 5
1. Compile the input generation code as follows : ```g++ -fopenmp parallel_5.cpp```

2. To generate 1000 bodies, using 8 chunks, and if you want to
   parallelize using 32 threads, please run : ```./a.out 1000.txt 8 32```
   
3. This will execute the code using a team of 8 threads which again create 32 threads
   for parallelization.