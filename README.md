# Heninger-Shacham Algorithm Implementation

We aimed to implement the Heninger-Shacham algorithm as detailed in ["Reconstructing RSA Private Keys from Random Key Bits"](https://link.springer.com/chapter/10.1007/978-3-642-03356-8_1). This algorithm is partially based on the premise that an attacker can exploit a physical property of DRAM memories by conducting a [Cold Boot Attack](https://citp.princeton.edu/our-work/memory/) and performing a memory dump.

DRAMs do not erase bits stored within them immediately after shutdown. Instead, the capacitors leak charges over time until reaching their ground states, opening a wide window for exploitation.

## Components

The implementation comprises mainly two parts:

### Heninger_part1.c

- **Functionality:** Deals with finding the value of K, correcting the most significant half of d, and finally computing kp & kq. We used the [Tonelli-Shanks algorithm](https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm) to compute the roots of the equation.
- **Details:** For exhaustive details, please refer to the mentioned paper.

### Heninger_part2.c

- **Functionality:** Handles the branch and prune part, verifying all possible combinations of bits for the components (d, p, q, d_p, d_q) at a given position "bit i". If they verify certain conditions mentioned in the paper, they are added in the position i; otherwise, they are pruned.
- **Note:** The keys must be in a specific format for `Heninger_part2.c` to work properly.

## Compilation and Execution

### Generating Keys

1. Compile `generate_key.c`:

```bash
gcc -o outputName generate_key.c -lgmp
```
