/**
 * @fileoverview Copyright (c) 2024-20xx, Stefano Gualandi,
 *               via Ferrata, 5, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@unipv.it (Stefano Gualandi)
 *
 */

#include "Utils.h"
#include "ATSP.h"

#include "HardTSPLIB_n_7_19.h"

void random_atsp(int n_max = 20)
{
  vector<double> C;
  C.reserve(n_max);

  for (int n = 10; n < 1 + n_max; n++)
  {
    C.clear();
    double opt = 0.0;

    ATSP solver(n);

    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        if (i != j)
          C.push_back((int)rand() % 512);

    auto start = std::chrono::high_resolution_clock::now();
    fprintf(stdout, "\nSolving ATSP with %d nodes:\t", n);
    solver.solveTSP(C, &opt);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = getMs(start, end);
    fprintf(stdout, "TSP: n: %d, opt: %.3f, runtime: %.3f\n", n, opt, elapsed);
  }
}

void test_ftv38_atsp()
{
  int n = int(39);
  ATSP solver(n);
  double tsp;
  vector<double> sol = solver.solveTSP(ftv38_atsp, &tsp);
}

void test_br17_atsp()
{
  int n = int(17);
  ATSP solver(n);
  double tsp;
  vector<double> sol = solver.solveTSP(br17_atsp, &tsp);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    fprintf(stdout, "run random instances for n = 10, ..., 100\n");
    random_atsp(100);
  }
  else if (argc == 2)
  {
    // Two tests from TSPLIB
    test_br17_atsp();
    test_ftv38_atsp();

    // Small hard instances from Hard-TSPLIB
    int nn = 7;
    double tsp;
    for (vector<double> &atsp : small_hard_atsp)
    {
      int n = nn;
      nn++;
      ATSP solver(n);
      vector<double> sol = solver.solveTSP(atsp, &tsp);
    }
  }

  return EXIT_SUCCESS;
}