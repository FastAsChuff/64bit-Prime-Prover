# 64bit-Prime-Prover
Trial division based deterministic primality test of 64 bit integers.

Test if argument is definitely prime using naive trial division. For large arguments, the program may use trial division by primes if primegaps.dat is available. This file is generated from primes.txt (all 203280221 32bit primes) if present. See https://github.com/FastAsChuff/Primes-List to create primes.txt = all 32bit primes = output of primeslist.bin 203280221

Memory/Storage Requirement: 250MB (if files used), <1MB otherwise.

Usage:- ./isprime2.bin number

0 <= number < 2^64

