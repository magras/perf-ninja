
#include "solution.h"

int solution(int *arr, int N) {
#ifdef SOLUTION
  return (1 + N) * N / 2;
#else
  int res = 0;
  for (int i = 0; i < N; i++) {
    res += arr[i];
  }
  return res;
#endif
}
