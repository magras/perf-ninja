
#include <array>

// Assume those constants never change
constexpr int N = 10000;
constexpr int minRandom = 0;
constexpr int maxRandom = 100;

// FIXME: this data structure can be reduced in size
struct S {
#ifdef SOLUTION
  int i : 8;
  short s : 8;
  long long l : 15;
  bool b : 1;
  float d;
#else
  int i;
  long long l;
  short s;
  double d;
  bool b;
#endif

  bool operator<(const S &s) const { return this->i < s.i; }
};

void init(std::array<S, N> &arr);
S create_entry(int first_value, int second_value);
void solution(std::array<S, N> &arr);
