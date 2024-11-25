#define SMALLPRIMETEST 100

typedef struct {
  uint32_t primedivchecksize;
  uint32_t *p;
  uint64_t *a, *m;
} primedivcheck_t;


void primedivcheckfree(primedivcheck_t *primedivcheck) {
  free(primedivcheck->p);
  primedivcheck->p = NULL;
  free(primedivcheck->a);
  free(primedivcheck->m);
} 

typedef struct {
  uint64_t w,L,count,sum;
  uint64_t *array;
} isprime5exceptions_t;

_Bool loadexceptions(isprime5exceptions_t *ex, char *filename) {
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) return false;
  uint64_t prev, temp, count = 0, sum = 0;
  if (fscanf(fp, "%lu %lu %lu %lu", &ex->w, &ex->L, &ex->count, &ex->sum) != 4) {
    fclose(fp);
    return false;
  }
  prev = 0;
  if ((ex->count == 0) || (ex->L >= 32) || (ex->w < 2)) {
    fclose(fp);
    return false;
  }
  ex->array = malloc(ex->count * sizeof(uint64_t));
  for (uint32_t i=0; i<ex->count; i++) {
    if (fscanf(fp, "%lu", &ex->array[i]) != 1) {
      fclose(fp);
      free(ex->array);
      return false;
    }
    if (ex->array[i] <= prev) {
      fclose(fp);
      free(ex->array);
      return false;
    }
    count++;
    sum += ex->array[i];
    prev = ex->array[i];
  }
  if ((sum != ex->sum) || (fscanf(fp, "%lu", &temp) == 1)) {
    fclose(fp);
    free(ex->array);
    return false;
  }
  fclose(fp);
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

_Bool inascarrayu64(uint64_t n, uint64_t *array, uint64_t arraysize) {
  // Binary Search. Array is assumed sorted ascending.
  if ((n < array[0]) || (n > array[arraysize-1])) return false;
  uint64_t leftix = 0;
  uint64_t rightix = arraysize-1;
  if (array[leftix] == n) return true;
  if (array[rightix] == n) return true;
  uint64_t middleix = (leftix + rightix)/2;
  while (true) {
    if (array[middleix] == n) return true;
    if (array[middleix] < n) {
      leftix = middleix;
    } else {
      rightix = middleix;
    }
    middleix = (leftix + rightix)/2;
    if (middleix == leftix) return false;
  }  
  return false;
}

_Bool isprime64dce(uint64_t n, primedivcheck_t primedivcheck, isprime5exceptions_t ex) {
  // Deterministic primality test using divisor check and exceptions list.
  // =====================================================================
  // Exceptions list ex.array[] is ALL composite 64 bit unsigned integers N such that
  //   N has no prime factors less than 2^(ex.L)
  // and
  //   (ex.w)^((N-1)/2) = +- 1 mod N
  // Therefore, given arbitrary odd 64 bit unsigned integer N > 2, 
  //   N is composite iff at least one of the following is true,
  //     (ex.w)^((N-1)/2) != +- 1 mod N
  //     N has a prime factor less than min(2^(ex.L), 1+isqrt(N)), 
  //     N is in ex.array[]
  if (n<=1) return false;
  if (n <= 3) return true;
  if ((n & 1) == 0) return false;
  uint64_t pmax = isqrt(n);
  if (n == pmax*pmax) return false; 
  uint32_t i = 0;
  for (; ((i < SMALLPRIMETEST) && (i < primedivcheck.primedivchecksize)); i++) {
    if (primedivcheck.p[i] > pmax) return true;
    if (primedivcheck.a[i]*n <= primedivcheck.m[i]) return false;
  }
  uint64_t res = modpowu64(ex.w, (n-1)/2, n);
  if ((res != 1) && (res != n-1)) return false;
  for (; (i < primedivcheck.primedivchecksize); i++) {
    if (primedivcheck.p[i] > pmax) return true;
    if (primedivcheck.a[i]*n <= primedivcheck.m[i]) return false;
  }
  // Better to use Hash Table instead for small ex.L.
  if (inascarrayu64(n, ex.array, ex.count)) return false;
  return true;
}

_Bool nextprime64dce(uint64_t n, uint64_t *nextprime, primedivcheck_t primedivcheck, isprime5exceptions_t ex) {
  while (n < 0xffffffffffffffffULL) {
    n++;
    if (isprime64dce(n, primedivcheck, ex)) {
      *nextprime = n;
      return true;
    }
  }
  return false;
}

_Bool loadprimes(isprime5exceptions_t ex, primedivcheck_t *primedivcheck, char *filename) {  
  primedivcheck->p = NULL;
  if (ex.array == NULL) return false;
  uint32_t primeslt2ton[] = {0, 0, 2, 4, 6, 11, 18, 31, 54, 97, 172, 309, 564, 1028, 1900, 3512, 
    6542, 12251, 23000, 43390, 82025, 155611, 295947, 564163, 1077871, 
    2063689, 3957809, 7603553, 14630843, 28192750, 54400028, 105097565, 203280221};
  uint64_t sumprimeslt2ton[] = {0, 0, 5, 17, 41, 160, 501, 1720, 6081, 22548, 80189, 289176, 1070091, 3908641, 14584641, 54056763, 202288087, 761593692, 2867816043, 10862883985, 41162256126, 156592635694, 596946687124, 2280311678414, 8729068693022, 33483086021512, 128615914639624, 494848669845962, 1906816620981654, 7357074544482779, 28422918403819825, 109930816131860852, 425649736193687430};
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) return false;
  uint64_t primesum = 2;
  uint32_t prime;
  uint32_t prevprime = 2;
  if (fscanf(fp, "%u", &prime) != 1) {
    fclose(fp);
    return false;
  }
  if (prime != 2) {
    fclose(fp);
    return false;
  }
  primedivcheck->primedivchecksize = primeslt2ton[ex.L]-1; // Odd primes only
  primedivcheck->p = malloc(primedivcheck->primedivchecksize*sizeof(uint32_t));
  for (uint32_t i=1; i<primeslt2ton[ex.L]; i++) {
    if (fscanf(fp, "%u", &prime) != 1) {
      free(primedivcheck->p);
      primedivcheck->p = NULL;
      fclose(fp);
      return false;
    }
    if (prime <= prevprime) {
      free(primedivcheck->p);
      primedivcheck->p = NULL;
      fclose(fp);
      return false;
    }
    primesum += prime;
    primedivcheck->p[i-1] = prime;
    prevprime = prime;
  }
  if (primesum != sumprimeslt2ton[ex.L]) {
    free(primedivcheck->p);
    primedivcheck->p = NULL;
    fclose(fp);
    return false;
  }
  primedivcheck->a = malloc(primedivcheck->primedivchecksize*sizeof(uint64_t));
  if (primedivcheck->a == NULL) {
    free(primedivcheck->p);
    primedivcheck->p = NULL;
    fclose(fp);
    return false;
  }
  primedivcheck->m = malloc(primedivcheck->primedivchecksize*sizeof(uint64_t));
  if (primedivcheck->m == NULL) {
    free(primedivcheck->p);
    primedivcheck->p = NULL;
    free(primedivcheck->a);
    fclose(fp);
    return false;
  }
  for (uint32_t i=0; i<primedivcheck->primedivchecksize; i++) {
    prime = primedivcheck->p[i];
    primedivcheck->a[i] = modinv64x(prime);
    primedivcheck->m[i] = 0xffffffffffffffffULL / prime;
  }
  return true;
}

