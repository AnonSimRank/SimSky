# SimSky

---

## Instruction of Matlab implementation

SimSky.m : full version, computing single-source SimRank search, includes the approximation of the diagonal correction matrix.   
SimSky_.m : simplified version, the default diagonal correction matrix is the identity matrix. 


**Matlab version:** Matlab R2020a

### Download code

```
git clone https://github.com/AnonSimRank/SimSky.git
cd SimSky/mcode
```

### Example

**Example1.m:** Test the ApproDiag algorithm on a graph containing six nodes.

**Example2.m:** Test the SimSky_ algorithm.

**Example3.m:** Test the actual error and upper bound.

**Example4.m:** Test the SimSky algorithm.

---

## Instruction of c++ implementation

**Compile Environment:** Visual Studio 2019 + MSVC16

**Dependence:**

   1. Eigen 

### Download code

```
git clone https://github.com/AnonSimRank/SimSky.git
cd SimSky/cppcode
```

### Example

The test data is included in file test.txt, put it in the same folder of SimRank.exe and run this program.

#### Compile and Run

1. Compile the project with MSVC16, select release|x64 configuration.  
2. Run the compiled program.  
