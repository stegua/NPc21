/**
 * @fileoverview Copyright (c) 2024-20xx, Stefano Gualandi,
 *               via Ferrata, 5, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@unipv.it (Stefano Gualandi)
 *
 */
#pragma once

#include <vector>
using std::vector;

#include <algorithm>
#include <chrono>

extern "C"
{
#include "gurobi_c.h"
#include "math.h"
}

#define POST(s, x)                                                 \
  if (x)                                                           \
  {                                                                \
    fprintf(stdout, "KA KA KA BOOM!!! %s: Error code %d\n", s, x); \
    exit(1);                                                       \
  }

const double INT_TOL = 0.00001;
const double SCALE = 1e+06;

// trasform time points diff to milliseconds
template <typename T>
inline double getMs(T start, T end)
{
  return double(
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                 .count()) /
         1000;
}

// Read a Text file and store it in a unique string
std::string readText(const std::string &filename, size_t &len)
{
  // Adapted from: https://github.com/simdjson/simdjson
  std::FILE *fp = std::fopen(filename.data(), "rb");

  if (fp == nullptr)
    return "";

  // Get the file size
  if (std::fseek(fp, 0, SEEK_END) < 0)
  {
    std::fclose(fp);
    return "";
  }

  long llen = std::ftell(fp);
  if ((llen < 0) || (llen == LONG_MAX))
  {
    std::fclose(fp);
    return "";
  }

  // Allocate the string
  len = static_cast<size_t>(llen);
  std::string s;
  s.reserve(len);
  if (s.data() == nullptr)
  {
    std::fclose(fp);
    return "";
  }

  // Read the string
  std::rewind(fp);
  size_t bytes_read = std::fread((char *)(s.data()), 1, len, fp);

  // TODO: strange character at the end
  // Found a solution at
  // https://stackoverflow.com/questions/17707904/fread-storing-random-characters-in-buffer
  // It seems that we need to add a NULL character at the end of the reading file
  // s[bytes_read] = '\0';
  if (std::fclose(fp) != 0 || bytes_read != len)
    return "";

  return s;
}
