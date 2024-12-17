// modpowu64() here https://github.com/FastAsChuff/Fast-Modular-Exponentiation
//=============================================================================================================
// Author: Simon Goater Dec 2024
int jacobiu64ss(unsigned long long a, unsigned long long n) { 
  // a,n > 0, n odd and gcd(a,n) = 1
  a %= n;
  if (n == 1ULL) return 1;
  if (a == 1ULL) return 1;
  if (a == 2) return 1 - 2*(int)(1ULL & ((n*n - 1ULL) >> 3));
  if (~a & 1ULL) return jacobiu64ss(2,n)*jacobiu64ss(a/2,n);
  return (1 - 2*(int)(1ULL & ((a - 1ULL)*(n - 1ULL) >> 2)))*jacobiu64ss(n, a);
}

_Bool isprimeu64ss(unsigned long long n) {
  // Solovay-Strassen primality test.
  // Deterministic for u64 inputs provided Jan Feitsma's file of u64 base 2 pseudoprimes is complete and correct.
  if (n < 2) return false;
  unsigned int bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
  unsigned int i, numbases = sizeof(bases)/sizeof(bases[0]);
  for (i=0; i<numbases; i++) {
    if (n == bases[i]) return true;
    if ((n % bases[i]) == 0) return false;
  }
  unsigned long long pow;
  int b;
  for (i=0; i<numbases; i++) {
    pow = modpowu64(bases[i], (n-1)/2, n);
    if ((pow != 1) && (pow != n-1)) return false;
    b = 1 - 2*(pow != 1);
    if (jacobiu64ss(bases[i], n) != b) return false;
  }
  unsigned long long exceptions[] = {
    91230634325542321ULL, 129545102216217601ULL, 289980482095624321ULL, 525025434548260801ULL,
    1606205228509922041ULL, 2238595159225471201ULL, 5322181695098476321ULL, 5830537550935874401ULL,
    6730931980364407201ULL, 6840760928649624001ULL, 9284340229544671201ULL, 14124830371276779001ULL};
  unsigned int numexceptions = sizeof(exceptions)/sizeof(exceptions[0]);
  for (i=0; i<numexceptions; i++) {
    if (exceptions[i] == n) return false;
  }
  return true;
}
//=============================================================================================================
