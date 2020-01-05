# Rust library for the course "Finite-volume methods" 

## Completed assignments 
* **a1** --- Godunov's method for linear PDEs (1D) 
* **a2** --- Exact Riemann solver for the Euler equations (1D) 
* **a3** --- Approximate Riemann solvers (1D) 
* **a4** --- Second-order methods with linear reconstruction and predictor-corrector time marching.   
* **a5-part1** --- First order 2D Euler equations solver using the HLLE approximate flux function.  
* **a5-part2** --- Second order 2D Euler equations solver using piecewise linear reconstruction and predictor-corrector time marching.   

Build all assignments using `cargo build --release`

Run assignments one-by-one using 
```shell
cargo run --bin <assignment_name> --release
# e.g. for assignment 1, 
cargo run --bin a1 --release
# assignment 5 is split between two executables: 
cargo run --bin a5-part1 --release
cargo run --bin a5-part2 --release
```

Run the tests with `cargo test`
