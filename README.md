Implementation of zkNARK for R1CS and its split accumulation scheme,
as described in [Proof-Carrying Data without Succinct Arguments](https://eprint.iacr.org/2020/1618.pdf).

## Quick Benchmarks

### zkNARK for R1CS
On a M1 Pro MacBook Pro

| Number of constraints | Proving time |
| --- | --- |
| 2^10 | 177.84 ms  |
| **2^13** | **1.2329 s** |
| 2^15 | 4.8807 s |

The circuit used in spartan-ecdsa consists of 8,076 constraints ~= 2^13.

## Run tests
```
cargo test
```

## Run benchmarks
```
cargo bench
```