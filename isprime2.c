#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

// gcc isprime2.c -o isprime2.bin -lm -O3 -march=native -Wall
#define NUMOF32BITPRIMES 203280221
char* primesfilename = "primes.txt";
char* primegapsfilename = "primegaps.dat";
// Bytes are gap/2 - 1 for gaps after 3. NUMOF32BITPRIMES - 2 bytes.
// sha256(primegapsfilename) = 2f38d96faadb2802723d90e4d9d240bef35a0f5092430794804fd9a75bbf25ed

_Bool loadprimegaps(uint8_t **primegaps, char *filename) {
  FILE *fp = fopen(filename, "rb");
  if (fp == NULL) return false;
  *primegaps = malloc(NUMOF32BITPRIMES-2);
  if (*primegaps == NULL) return false;
  size_t bytes = fread(*primegaps, 1, NUMOF32BITPRIMES-2, fp);
  fclose(fp);
  if (bytes != NUMOF32BITPRIMES-2) return false;
  return true;
}

_Bool saveprimegaps(uint8_t *primegaps, char *filename) {
  FILE *fp = fopen(filename, "wb");
  if (fp == NULL) return false;
  size_t bytes = fwrite(primegaps, 1, NUMOF32BITPRIMES-2, fp);
  fclose(fp);
  if (bytes != NUMOF32BITPRIMES-2) return false;
  return true;
}

uint8_t *loadprimegaps0(char *filename) {
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) return NULL;
  uint64_t primesum = 5;
  uint32_t prevprime = 3;
  uint32_t prime, diff;
  if (fscanf(fp, "%u", &prime) != 1) return NULL;
  if (prime != 2) return NULL;
  if (fscanf(fp, "%u", &prime) != 1) return NULL;
  if (prime != 3) return NULL;
  uint8_t *res = malloc(NUMOF32BITPRIMES-2);
  if (res == NULL) return NULL;
  uint32_t i = 0;
  while (fscanf(fp, "%u", &prime) == 1) {
    diff = prime - prevprime;
    if ((diff > 512) || (i >= NUMOF32BITPRIMES-2)) {
      free(res);
      fclose(fp);
      return NULL;
    }
    res[i] = (diff/2) - 1;
    prevprime = prime;
    primesum += prime;
    i++;
  }
  fclose(fp);
  //printf("prime sum is %lu\n", primesum);
  if (primesum != 425649736193687430ULL) {
    free(res);
    return NULL;
  }
  return res;
}

uint32_t isqrt(uint64_t n) {
  if (n < 2) return n;
  uint64_t ai = sqrt(n);
  while (!((ai <= n/ai) && ((ai+1) > n/(ai+1)))) {    
    ai = (ai + n/ai)/2;
  }
  return ai;
}

_Bool isprime64g(uint64_t n, uint8_t *primegaps) {
  // Deterministic primality test using trial division by primes.
  if (n<=1) return false;
  if (n <= 3) return true;
  if (n == 5) return true;
  if ((n & 1) == 0) return false;
  if ((n % 3) == 0) return false; 
  if ((n % 5) == 0) return false; 
  uint64_t pmax = isqrt(n);
  if (n == pmax*pmax) return false; 
  uint32_t prime = 5;
  for (uint32_t i = 1; i < NUMOF32BITPRIMES-2; i++) {
    prime += (((uint32_t)primegaps[i] + 1U)*2);
    if (prime > pmax) break;
    if ((n % prime) == 0) return false;
  }
  return true;
}

_Bool nextprime64g(uint64_t n, uint64_t *nextprime, uint8_t *primegaps) {
  while (n < 0xffffffffffffffffULL) {
    n++;
    if (isprime64g(n, primegaps)) {
      *nextprime = n;
      return true;
    }
  }
  return false;
}

_Bool isprime64(uint64_t n) {
  // Deterministic primality test using trial division.
  if (n<=1) return false;
  if (n <= 3) return true;
  if (n == 5) return true;
  if ((n & 1) == 0) return false;
  if ((n % 3) == 0) return false;
  if ((n % 5) == 0) return false;
  uint32_t imax = isqrt(n);
  // 1,7,11,13,17,19,23,29 mod 30 only
  for (uint64_t i = 7; ; i+=6) {
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 4;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 2;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 4;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 2;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 4;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 6;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 2;
    if (n%i == 0) return false;
    if (i >= imax) break;
  }
  return true;
}

_Bool nextprime64(uint64_t n, uint64_t *nextprime) {
  while (n < 0xffffffffffffffffULL) {
    n++;
    if (isprime64(n)) {
      *nextprime = n;
      return true;
    }
  }
  return false;
}

uint64_t atou64(char *in) {
  uint64_t res = 0;
  while (*in) {
    res *= 10;
    res += *in - '0';
    in++;
  }
  return res;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Test if argument is definitely prime using naive trial division. For large arguments, the program may use trial division by primes if %s is available. This file is generated from %s (all 203280221 32bit primes) if present.\nMemory/Storage Requirement: 250MB (if files used), <1MB otherwise.\nUsage:- %s number\n0 <= number < 2^64\n", primegapsfilename, primesfilename, argv[0]);
    exit(1);
  }
  uint8_t *primegaps = NULL;
  uint64_t n = atou64(argv[1]);
  if (n > 1000000000000000ULL) {
    if (!loadprimegaps(&primegaps, primegapsfilename)) {
      printf("Load of %s failed.\n", primegapsfilename);
      // See https://github.com/FastAsChuff/Primes-List
      // Create primes.txt = all 32bit primes = output of primeslist.bin 203280221
      primegaps = loadprimegaps0(primesfilename);
      if (primegaps) {
        if (!saveprimegaps(primegaps, primegapsfilename)) {
          printf("Save of %s failed.\n", primesfilename);          
        }
      } else {
        printf("Load of %s failed.\n", primesfilename);  
      }
    }
  }
  uint64_t next;
  _Bool nisprime;
  if (primegaps) {
    nisprime = isprime64g(n, primegaps);
  } else {
    nisprime = isprime64(n);
  }
  if (nisprime) {
    printf("%lu is proven prime.\n", n);
  } else {
    printf("%lu is not prime.\n", n);
  }
  if (primegaps) {
    if (nextprime64g(n, &next, primegaps)) printf("Next proven prime is %lu\n", next);
  } else {
    if (nextprime64(n, &next)) printf("Next proven prime is %lu\n", next);
  }
}
