# 64bit-Prime-Prover
Trial division based deterministic primality test of 64 bit integers.

Test if argument is definitely prime using naive trial division. For large arguments, the program may use trial division by primes if primegaps.dat is available. This file is generated from primes.txt (all 203280221 32bit primes) if present. See https://github.com/FastAsChuff/Primes-List to create primes.txt = all 32bit primes = output of primeslist.bin 203280221

The main motive for writing these programs was out of curiosity about just how fast, or slow, trial division can be on 64 bit inputs. It was not intended to be practical, but may be useful to anyone who doesn't want to rely on deep number theory, opaque library functions, or claims of exhaustive testing which are extrememly difficult to verify.

isprime2.bin
------------

Memory/Storage Requirement: 250MB (if files used), <1MB otherwise.

Usage:- ./isprime2.bin number

0 <= number < 2^64

e.g. ./isprime2.bin 18361375334787048247

18361375334787048247 is proven prime.

Next proven prime is 18361375334787048287

isprime4.bin
------------

Memory/Storage Requirement: 4GB (if files used), <1MB otherwise.

Usage:- ./isprime4.bin number

0 <= number < 2^64

e.g. ./isprime4.bin 18446744073680551637

Time of test = 625 ms

Next proven prime is 18446744073680551639

