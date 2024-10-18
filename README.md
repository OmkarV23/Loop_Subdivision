# Loop Subdivision

This project is an implementation of the Loop Subdivision algorithm for triangular meshes. The algorithm is based on the paper "Smooth Subdivision Surfaces Based on Triangles" by Charles Loop. The implementation is done in C++.

## Build
To build the project, run the following commands:
```bash
cd Subdivision
make
```

## Run
To run the project, execute the following command:
```bash
./Subdivision <input_file> <output_file>
```

## Stanford Bunny
<p align="center">
  <img src="/Results/bunny_base.png" alt="Iteration 0" width="700">
  <!-- <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/32/White_Right_Arrow.svg/120px-White_Right_Arrow.svg.png" alt="Arrow" width="50" height="50"> -->
  <img src="/Results/bunny_lsd_1.png" alt="Iteration 1" width="700">
  <!-- <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/32/White_Right_Arrow.svg/120px-White_Right_Arrow.svg.png" alt="Arrow" width="50" height="50"> -->
  <img src="/Results/bunny_lsd_2.png" alt="Iteration 2" width="700">
</p>

## Utah Teapot
<p align="center">
  <img src="/Results/teapot_base04.png" alt="Iteration 0" width="700">
  <!-- <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/32/White_Right_Arrow.svg/120px-White_Right_Arrow.svg.png" alt="Arrow" width="50" height="50"> -->
  <img src="/Results/teapot_lsd_1.png" alt="Iteration 1" width="700">
  <!-- <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/32/White_Right_Arrow.svg/120px-White_Right_Arrow.svg.png" alt="Arrow" width="50" height="50"> -->
  <img src="/Results/teapot_lsd_2.png" alt="Iteration 2" width="700">
</p>