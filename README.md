
## Aim
A simple program to convert an input Nastran punch and an analysis input files into ASCII files

## Compilation
The program can be easily compiled in the following way: g++ -o PchConverter PchConverter.cpp -O3 -std=c++17

## Input files
- Nastran punch file (with '.pch')
- Analysis input file (with '.dat')

## Result files
- Sparse stiffness matrix
- Sparse mass matrix
- Sparse damping matrix
- Mapping file
- Nodes coordinates

## User interface
There are two ways to specify the input files:
1. The file loader with the name "Loader.txt" ($FILE_LOADER_NAME) located in the directory with the executable. Its content has to be formatted in the following way:
```
<Nastran punch file>
<Analysis input file>
<Geometry scale> 
```
For example,
```
examples/beam.pch
examples/beam.dat
1.0
```
2. User input