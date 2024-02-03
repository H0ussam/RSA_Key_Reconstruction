We aimed to implement the Heninger-shacham algorithm [Reconstructing RSA Private Keys from Random Key Bits](https://link.springer.com/chapter/10.1007/978-3-642-03356-8_1), it is partially based on the premise that an attacker can exploit a phyicaly property of memories "type DRAM" by conducting a [Cold Boot Attack](https://citp.princeton.edu/our-work/memory/) and reforming a dump of the memory.
It turned out that DRAMs don't earease bits stored within it immediatly after shutdown, instead the capacitors leak charges over time until reaching thier ground states, this open a wide window for exploitation.

you find mainly 2 parts of the whole algorithm:

the file named Heninger_part1.c, deals with finding the value of K, correcting the most significant half of d, finally computing kp & kq, we used the [Tonelli-Shanks algorithm](https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm), to compute the roots of the equation.
For exhaustive details please read the mentioned paper above.

The file named  Heninger_part2.c deals witht he branch and prune part, it verifies all possible combination of bits for the coponents (d, p, q, d_p, d_q) at a given position "bit i", if they verifiy certain conditions mentioned in the paper, they are added to in the position i otherwise pruned; again for exhaustive details check the paper x)

The keys degraded should be in a certain format for Heninger_part2.c to ork properly

Compile the file generate_key.c:

'$gcc -o outputName generate_key.c -lgmp'

