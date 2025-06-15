# 1D-FFT-in-RISCV
1D FFT algorithm written in RISC V assembly through vectorization.
This project was created for Computer Architecture and Assembly Language course. The simulator that this code was run on was whisper and on a virtual machine (credits to @awsay905), hence the json and makefile.
Note that it is necessary to properly set the value of "bytes_per_vec" and "max_bytes_per_elem" to set the max number of elements in a vector, otherwise the algorithm will not run completely.
The functions for bit reversal were contributed by @akx16
