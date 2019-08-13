# Adyar-RS

Approximate Common Substring

## Dependencies

* [Rust](https://www.rust-lang.org/) compiler version >= 1.36. See install instructions [here](https://www.rust-lang.org/learn/get-started).
* A modern, C++11 ready compiler such as `g++` >= v4.7  or `clang` >= v3.2 or MSVC >= 14 to compile libdivsufsort.
* The [cmake](www.cmake.org) build system (Version >= 2.8.11).
* A 64-bit operating system. Tested in Windows 10 and Linux.

## Compiling the code

The program can be compiled as follows:

     cargo build

To create a release build

     cargo build --release

In case of release build, the output executable will be available at `target/release/adyar-rs`

## Running the program

The release version of he program is run one of the following two ways :
`cargo run --release -- <FASTA_FILE> -o <DIST_MATRIX> -k <K_VALUE> `

or

`target/release/adyar-rs <FASTA_FILE> -o <DIST_MATRIX> -k <K_VALUE> `

Example:
`target/release/adyar-rs data/input.fa -o input_dist.txt -k 2`

Output:
Runtime, value of symmetric acs distance
