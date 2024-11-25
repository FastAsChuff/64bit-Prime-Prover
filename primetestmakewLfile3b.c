#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include </home/simon/fastmodinvpow2fns.c> // See https://github.com/FastAsChuff/Fast-Modular-Inverse-Modulo-Powers-Of-2
#include </home/simon/modpowu64.c> // See https://github.com/FastAsChuff/Fast-Modular-Exponentiation/tree/main

// gcc primetestmakewLfile3b.c -o primetestmakewLfile3b.bin -lm -pthread -O3 -march=native -Wall

#define NUMOF32BITPRIMES 203280221
#define NUMOF16BITPRIMES 6542
#define MAXTHREADS 1000
#define PRIMESPERTHREAD 2500
#define SMALLPRIMETEST 50

void die(char *str) {
  fprintf(stderr, "%s\n", str);
  exit(1);
}

uint32_t modpowu32(uint32_t a, uint32_t e, uint32_t n) {
// Returns a^e mod n
  if (n == 0) return 0;
  if (a < 2) return a;
  uint64_t res = 1;
  uint64_t sq = a % n;
  while (e) {
    if (e & 1U) res = (res * sq) % n;
    sq = (sq*sq) % n;
    e >>= 1;
  }
  return res;
}

typedef struct {
  uint32_t *p;
  uint64_t *a, *m;
} primedivcheck_t;

uint32_t multorderu32(uint32_t a, uint32_t p, primedivcheck_t primedivcheck) {
// For prime p (primality assumed), returns least +ve m s.t. a^m = 1 mod p
// [Compared results with naive function for 2 <= a < 1000 and all 16 bit primes.]
  if (p < 2) return 0;
  a = a % p;
  if (a < 2) return a;
  uint32_t e = p-1;
  uint32_t primefactors[12];
  uint32_t primefactorsfound = 0;
  uint32_t inc = 0;
  while ((e & 1) == 0) {
    inc = 1;
    primefactors[0] = 2;
    e /= 2;
  }
  primefactorsfound += inc;
  uint32_t i = 0;
  while ((i < (NUMOF16BITPRIMES-1)) && (e > 1)) {
    inc = 0;
    while (e*primedivcheck.a[i] <= primedivcheck.m[i]) {
      inc = 1;
      e /= primedivcheck.p[i];
      primefactors[primefactorsfound] = primedivcheck.p[i];
    }
    primefactorsfound += inc;
    i++;
  }
  if (e > 1) primefactors[primefactorsfound++] = e;
  e = p-1;
  for (i = 0; i<primefactorsfound; i++) {
    while (modpowu32(a,e/primefactors[i],p) == 1) {
      e /= primefactors[i];
      if (e % primefactors[i] != 0) break;
    }
  }
  return e;
}

char* primesfilename = "primes.txt";

uint64_t atou64(char *in) {
  uint64_t res = 0;
  while (*in) {
    res *= 10;
    res += *in - '0';
    in++;
  }
  return res;
}

_Bool load32bitprimes(primedivcheck_t *primedivcheck, uint64_t twotoL, uint32_t *primedivchecksize, char *filename) {
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) return false;
  uint64_t primesum = 2;
  if (twotoL > 0x100000000ULL) return false;
  uint32_t prime;
  if (fscanf(fp, "%u", &prime) != 1) {
    fclose(fp);
    return false;
  }
  if (prime != 2) {
    fclose(fp);
    return false;
  }
  primedivcheck->p = malloc((NUMOF32BITPRIMES-1)*sizeof(uint32_t));
  if (primedivcheck->p == NULL) {
    fclose(fp);
    return false;
  }
  uint32_t i = 0;
  *primedivchecksize = 0;
  while (fscanf(fp, "%u", &prime) == 1) {
    if (i < NUMOF32BITPRIMES-1) primedivcheck->p[i] = prime;
    if (prime < twotoL) (*primedivchecksize)++;
    primesum += prime;
    i++;
  }
  fclose(fp);
  if (primesum != 425649736193687430ULL) {
    free(primedivcheck->p);
    primedivcheck->p = NULL;
    return false;
  }
  primedivcheck->a = malloc(*primedivchecksize*sizeof(uint64_t));
  if (primedivcheck->a == NULL) {
    free(primedivcheck->p);
    primedivcheck->p = NULL;
    return false;
  }
  primedivcheck->m = malloc(*primedivchecksize*sizeof(uint64_t));
  if (primedivcheck->m == NULL) {
    free(primedivcheck->p);
    primedivcheck->p = NULL;
    free(primedivcheck->a);
    return false;
  }
  for (i = 0; i<*primedivchecksize; i++) {
    primedivcheck->a[i] = modinv64x(primedivcheck->p[i]);
    primedivcheck->m[i] = 0xffffffffffffffffULL / primedivcheck->p[i];
  }
  return true;
}

uint32_t isqrt(uint64_t n) {
  if (n < 2) return n;
  uint64_t ai = sqrt(n);
  while (!((ai <= n/ai) && ((ai+1) > n/(ai+1)))) {    
    ai = (ai + n/ai)/2;
  }
  return ai;
}

int uint64asc(const void*a, const void*b) {
  const uint64_t *x = a;
  const uint64_t *y = b;
  if (*x == *y) return 0;
  return 2*(*x > *y) - 1;
}

typedef struct {  
  uint32_t threadno;
} findexceptionsargs_t;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
uint64_t primeix = 0;
pthread_t *threads;
findexceptionsargs_t *fexargs;
uint64_t *count;
uint64_t *countmax;
uint64_t **exceptions;
primedivcheck_t primedivcheck;
uint32_t primedivchecksize;
uint64_t numthreads, w = 2, L = 3;

void *findexceptions(void *findexceptionsargs) {
/*
For 22 <= L < 32, (mc+1) is prime for exceptions. Otherwise, n = p(mc+1) has a prime factor < 2^L.
p(mc+1) <= (2^64)-1
(mc+1) <= floor(((2^64)-1)/p)
m <= floor((floor(((2^64)-1)/p) - 1)/c)
To avoid finding the same n twice, can assume
p <= mc+1
floor((p+c-2)/c) = ceil((p-1)/c) <= m
*/
  findexceptionsargs_t *args = findexceptionsargs;
  uint32_t ixstart = primedivchecksize + args->threadno*PRIMESPERTHREAD;
  _Bool done = false;
  while (!done) {
    for (uint32_t i = ixstart; (i < NUMOF32BITPRIMES-1) && (i < ixstart + PRIMESPERTHREAD); i++) {
      uint32_t c = multorderu32(w % primedivcheck.p[i], primedivcheck.p[i], primedivcheck); 
      uint64_t mmax = ((0xffffffffffffffffULL/primedivcheck.p[i]) - 1)/c;
      uint64_t m = (((uint64_t)primedivcheck.p[i]+c)-2)/c;  // W.L.O.G mc+1 >= P => m >= ceil((P-1)/c)
      uint64_t mcp1;
      for (; m <= mmax; m++) {
        mcp1 = m*c + 1;
        if (mcp1 & 1ULL) {
          _Bool divisorfound = false;
          for (uint32_t j = 0; j < SMALLPRIMETEST; j++) {
            if (primedivcheck.a[j]*mcp1 <= primedivcheck.m[j]) {
              divisorfound = true;
              break;
            }
          }
          if (!divisorfound) {
              uint64_t n = primedivcheck.p[i]*mcp1;
              uint64_t wtonm1o2modn = modpowu64(w, (n-1)/2, n);
              if ((wtonm1o2modn == 1) || (wtonm1o2modn == n-1)) {
                divisorfound = false;
                uint64_t pmax = isqrt(mcp1);
                for (uint32_t j = SMALLPRIMETEST; j < primedivchecksize; j++) {
                  if (primedivcheck.p[j] > pmax) break;
                  if (primedivcheck.a[j]*mcp1 <= primedivcheck.m[j]) {
                    divisorfound = true;
                    break;
                  }
                }
                if (!divisorfound) {
                  if (count[args->threadno] >= countmax[args->threadno]) {
                    countmax[args->threadno] *= 2;
                    exceptions[args->threadno] = realloc(exceptions[args->threadno], countmax[args->threadno]*sizeof(uint64_t));
                    if (exceptions[args->threadno] == NULL) die("Failed to allocate enough memory.");
                  }
                  exceptions[args->threadno][count[args->threadno]] = n;
                  count[args->threadno]++;
                }
              }
          }
        }
      }
    }
    pthread_mutex_lock(&mutex);
    if (primeix < NUMOF32BITPRIMES-1) {
      ixstart = primeix;
      printf("Thread no %u continuing from prime ix %lu. %lu remaining.\n", args->threadno, primeix + PRIMESPERTHREAD, NUMOF32BITPRIMES-1 - primeix);
      primeix += PRIMESPERTHREAD;
    } else done = true;
    pthread_mutex_unlock(&mutex);
  }
  return NULL;
}

int main(int argc, char** argv) {
  _Bool inputvalid = true;
  char *outfilename = "primetestwLfile3.txt";
  if (argc >= 4) {
    numthreads = atou64(argv[1]);
    w = atou64(argv[2]);
    L = atou64(argv[3]);
  } else inputvalid = false;
  uint64_t twotoL = 1ULL << 30;
  if ((L > 32) || (L < 16)) {
    inputvalid = false;
  } else {
    twotoL = 1ULL << L;
  }
  if ((w >= twotoL) || (w < 2)) inputvalid = false;
  uint64_t totcount = 0;
  if ((numthreads > MAXTHREADS) || (numthreads == 0)) {
    inputvalid = false;
  } else {
    for (uint32_t i=0; i<numthreads; i++) {
      threads = malloc(numthreads*sizeof(pthread_t));
      if (threads == NULL) die("Failed to malloc threads!");
      fexargs = malloc(numthreads*sizeof(findexceptionsargs_t));
      if (fexargs == NULL) die("Failed to malloc fexargs!");
      count = malloc(numthreads*sizeof(uint64_t));
      if (count == NULL) die("Failed to malloc count!");
      countmax = malloc(numthreads*sizeof(uint64_t));
      if (countmax == NULL) die("Failed to malloc countmax!");
      exceptions = malloc(numthreads*sizeof(uint64_t*));
      if (exceptions == NULL) die("Failed to malloc *exceptions!");
    }
    for (uint32_t i=0; i<numthreads; i++) {
      count[i] = 0;
      countmax[i] = 10000; // Initial buffer size
      exceptions[i] = malloc(countmax[i]*sizeof(uint64_t));
      if (exceptions[i] == NULL) die("Failed to allocate enough memory.");
    }
  }
  uint64_t sum = 0;
  FILE *fp = fopen(outfilename, "w");
  if (fp == NULL) {
    printf("Failed to open %s for writing.\n", outfilename);
    inputvalid = false;
  }
  if (!inputvalid) {
    printf("This program outputs to %s all 64 bit unsigned composite integers N with no prime factors below 2^L, for which w^(N-1)/2 = +-1 mod N. The first line of the output file contains four positive integers, w L count sum, followed by the count composites found whose sum mod 2^64 is sum. The composites are output in strictly increasing order. Requires file %s of all 32 bit primes.\n", outfilename, primesfilename);
    printf("Usage: %s mumthreads w L\n", argv[0]);
    printf("16 <= L <= 31\n");
    printf("2 <= w < 2^L\n");
    printf("1 <= numthreads <= %u\n", MAXTHREADS);
    exit(0);
  }
  if (!load32bitprimes(&primedivcheck, twotoL, &primedivchecksize, primesfilename)) {
    printf("Error loading from %s\n", primesfilename);
    exit(1);
  }
  printf("Starting from %u\n", primedivcheck.p[primedivchecksize]);
  primeix = primedivchecksize + numthreads*PRIMESPERTHREAD;
  for (uint32_t i=0; i<numthreads; i++) {
    fexargs[i].threadno = i;
    pthread_create(&threads[i], NULL, findexceptions, &fexargs[i]);
  }
  for (uint32_t i=0; i<numthreads; i++) {
    pthread_join(threads[i], NULL);
  }
  for (uint32_t i=0; i<numthreads; i++) {
    totcount += count[i];    
    for (uint32_t j=0; j<count[i]; j++) sum += exceptions[i][j];
  }
  uint64_t *allexceptions =  malloc(totcount*sizeof(uint64_t));
  if (allexceptions == NULL) die("Failed to allocate enough memory.");
  uint64_t allexceptionsix = 0;
  for (uint32_t i=0; i<numthreads; i++) {
    memcpy(allexceptions + allexceptionsix, exceptions[i], count[i] * sizeof(uint64_t));
    allexceptionsix += count[i];
  }
  fprintf(fp, "%lu %lu %lu %lu\n", w, L, totcount, sum);
  qsort(allexceptions, totcount, sizeof(uint64_t), &uint64asc);
  for (uint64_t i=0; i<totcount; i++) {
    fprintf(fp, "%lu\n", allexceptions[i]);
  }
  fclose(fp);
  free(allexceptions);
  free(primedivcheck.p);
  free(primedivcheck.a);
  free(primedivcheck.m);
  for (uint32_t i=0; i<numthreads; i++) free(exceptions[i]);
  free(exceptions);
  free(countmax);
  free(count);
  free(threads);
  free(fexargs);
}
