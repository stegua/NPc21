/**
 * @fileoverview Copyright (c) 2024-20xx, Stefano Gualandi,
 *               via Ferrata, 5, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@unipv.it (Stefano Gualandi)
 *
 */

#pragma once

#include "Utils.h"

#include <lemon/adaptors.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>
#include <lemon/hao_orlin.h>
#include <lemon/smart_graph.h>

#include <lemon/connectivity.h>
using namespace lemon;

// Graph instance
typedef double Value;
typedef SmartDigraph Digraph;

typedef SmartDigraph::Node Node;
typedef SmartDigraph::Arc Arc;
typedef SmartDigraph::ArcMap<Value> CapMap;
typedef SmartDigraph::NodeMap<int> CutMap;

// Callback for LP relaxation
struct callback_data
{
  int n;
  int m;
  vector<int> head;
  vector<int> tail;

  // Pre-initialized graph and associated data
  SmartDigraph G;
  vector<SmartDigraph::Node> nodes;
  vector<SmartDigraph::Arc> arcs;

  // For adding cuts
  vector<double> rmatval;
  vector<int> rmatind;
};

int __stdcall subtour_elim(GRBmodel *model, void *cbdata, int where,
                           void *usrdata)
{
  // Get reference to data allocated only once
  struct callback_data *mydata = (struct callback_data *)usrdata;
  int m = mydata->m;
  const vector<int> &head = mydata->head;
  const vector<int> &tail = mydata->tail;
  SmartDigraph &Gn = mydata->G;
  vector<SmartDigraph::Node> &N = mydata->nodes;
  vector<SmartDigraph::Arc> &E = mydata->arcs;

  vector<double> &rmatval = mydata->rmatval;
  vector<int> &rmatind = mydata->rmatind;

  double *xbar = NULL;
  int error = 0;
  int status = 0;

  if (where == GRB_CB_MIPNODE)
  {
    xbar = (double *)malloc(m * sizeof(double));
    if (xbar == NULL)
    {
      fprintf(stderr, "Out of memory\n");
      return -1;
    }
    POST("get status: ", GRBcbget(cbdata, where, GRB_CB_MIPNODE_STATUS, &status));
    if (status == GRB_OPTIMAL)
    {
      POST("GRBcbget: ", GRBcbget(cbdata, where, GRB_CB_MIPNODE_REL, xbar));

      // Create an ArcMap to mark active/inactive edges
      SmartDigraph::ArcMap<bool> active(Gn, true);

      // Scaling to integer weights (WARNING: numerical tolerances, to check!)
      CapMap xb(Gn);
      for (int i = 0; i < m; i++)
        if (xbar[i] < INT_TOL)
          active[E[i]] = false; // Filter out the arc
        else
        {
          if (xbar[i] > 1 - INT_TOL)
            xb[E[i]] = SCALE;
          else
            xb[E[i]] = int(round(SCALE * xbar[i]));
        }

      // Create a filtered graph view
      FilterArcs<SmartDigraph> Gt(Gn, active);
      // Run over the filtered graph
      HaoOrlin<FilterArcs<SmartDigraph>, CapMap> ni(Gt, xb);
      ni.run();

      if (double(ni.minCutValue()) * 100 < 0.99 * SCALE)
      {
        CutMap cut(Gn);
        ni.minCutMap(cut);

        rmatind.clear();
        for (int i = 0; i < m; ++i)
          if (cut[N[head[i]]] && cut[N[head[i]]] != cut[N[tail[i]]])
            rmatind.push_back(i);

        if (rmatind.size() > 2)
          POST("add cbcut", GRBcbcut(cbdata, (int)rmatind.size(), &rmatind[0],
                                     &rmatval[0], GRB_GREATER_EQUAL, 1.0));
      }
    }
  }

  if (where == GRB_CB_MIPSOL)
  {
    xbar = (double *)malloc(m * sizeof(double));
    if (xbar == NULL)
    {
      fprintf(stderr, "Out of memory\n");
      return -1;
    }

    POST("GRBcbget GRB_CB_MIPSOL:", GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, xbar));

    // Create an ArcMap to mark active/inactive edges
    SmartDigraph::ArcMap<bool> active(Gn, false);

    // The solution must be integer
    for (int i = 0; i < m; ++i)
      if (xbar[i] > 0.5)
        active[E[i]] = true;

    // Create a filtered graph view
    FilterArcs<SmartDigraph> Gt(Gn, active);

    CutMap cut(Gn);

    int nn = lemon::stronglyConnectedComponents(Gt, cut);

    if (nn > 1)
    {
      for (int k = 0; k < nn; k++)
      {
        rmatind.clear();
        for (int i = 0; i < m; ++i)
          if (cut[N[head[i]]] == k && cut[N[tail[i]]] != k)
            rmatind.push_back(i);

        if (rmatind.size() > 2)
          POST("add lazy", GRBcblazy(cbdata, (int)rmatind.size(), &rmatind[0],
                                     &rmatval[0], GRB_GREATER_EQUAL, 1.0));
      }
    }
  }

  if (xbar != NULL)
    free(xbar);

  return error;
}

// Subtour polytope
class ATSP
{
public:
  ATSP(int _n) : n(_n), m(n * (n - 1))
  {
    // Set up nodes
    N.reserve(n);
    for (int i = 0; i < n; i++)
      N.emplace_back(mydata.G.addNode());
    // Set up edges
    E.reserve(m);
    Head.reserve(m);
    Tail.reserve(m);
    mydata.rmatind.reserve(m);
    mydata.rmatval.reserve(m);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        if (i != j)
        {
          E.push_back(mydata.G.addArc(N[i], N[j]));
          Head.push_back(i);
          Tail.push_back(j);
          mydata.rmatval.push_back(1);
        }

    mydata.n = n;
    mydata.m = m;
    mydata.head = Head;
    mydata.tail = Tail;
    mydata.nodes = N;
    mydata.arcs = E;

    // Define problem data
    vector<int> cmatcnt(m, 2);
    vector<int> cmatbeg(m, 0);
    vector<int> cmatind(m * 2, 0);
    vector<double> cmatval(m * 2, 1);

    vector<double> xbar(m, 0);

    size_t idx = 0;
    size_t e = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        if (i != j)
        {
          cmatbeg[e++] = idx;
          cmatind[idx++] = i;
          cmatind[idx++] = n + j;
        }

    // Variables c_ij = dist(i,j)
    vector<double> obj(m, 1);

    // Variables 0 <= x_ij <= 1
    vector<double> lower(m, 0);
    vector<double> upper(m, 1);

    // Constraints "sum_ij x_ij == 1"
    vector<double> rhs(2 * n, 1);
    vector<char> sense(2 * n, 'E');
    vector<char> vtype(m, 'B');

    GRBloadenv(&env_master, NULL);
    GRBsetintparam(env_master, GRB_INT_PAR_OUTPUTFLAG, 0);
    GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
    GRBsetintparam(env_master, GRB_INT_PAR_LAZYCONSTRAINTS, 1); // for cblazy callback
    GRBsetintparam(env_master, GRB_INT_PAR_PRECRUSH, 1);        // for cbcut callback

    // REMARKS: on larger instances, you can play with fine-tuning the parameters
    //          Otherwise, use the default values
    GRBsetintparam(env_master, GRB_INT_PAR_THREADS, 1);
    // GRBsetintparam(env_master, GRB_INT_PAR_MIPFOCUS, 3);
    // GRBsetintparam(env_master, GRB_INT_PAR_NUMERICFOCUS, 3);
    // GRBsetdblparam(env_master, GRB_DBL_PAR_OPTIMALITYTOL, 1e-09);

    fprintf(stdout, "n: %d, m: %d\n", n, m);
    POST("load model",
         GRBloadmodel(env_master, &master, "grb_atsp", m, 2 * n, 1, 0, &obj[0],
                      &sense[0], &rhs[0], &cmatbeg[0], &cmatcnt[0], &cmatind[0],
                      &cmatval[0], &lower[0], &upper[0], &vtype[0], NULL, NULL));
  }

  // Destructor
  ~ATSP()
  {
    if (master != nullptr)
      POST("free ATSP", GRBfreemodel(master));
    if (env_master != nullptr)
      GRBfreeenv(env_master);
  }

  // Solve the IP problem (i.e., the ATSP)
  vector<double> solveTSP(vector<double> &xb, double *optvalue)
  {
    // Change cost to the objective function
    POST("solveTSP::change cost coef:",
         GRBsetdblattrarray(master, "Obj", 0, (int)xb.size(), &xb[0]));

    int status;
    double value;
    double runtime;
    double nnnodes = 0;
    vector<double> xbar(m, 0);
    vector<Node> Nh;
    Nh.reserve(n);

    POST("add callback",
         GRBsetcallbackfunc(master, subtour_elim, (void *)&mydata));

    GRBoptimize(master);

    POST("get status: ", GRBgetintattr(master, "Status", &status));
    if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE)
    {
      fprintf(stdout, "Unbounded or Infeasible\n");
      exit(0);
    }

    POST("get X", GRBgetdblattrarray(master, "X", 0, m, &xbar[0]));
    POST("get obj", GRBgetdblattr(master, "ObjVal", &value));

    POST("get numnodes", GRBgetdblattr(master, "NodeCount", &nnnodes));
    POST("get runtime", GRBgetdblattr(master, "Runtime", &runtime));
    fprintf(stdout, "TSP status: %d => n: %d, fobj: %.2f, numnodes: %.f, runtime: %.2f\n",
            status, n, value, nnnodes, runtime);

    optvalue[0] = value;

    return xbar;
  }

private:
  int n;
  int m;

  // Pointers for Gurobi
  GRBenv *env_master;
  GRBmodel *master;

  // Graph data structure
  SmartDigraph G;
  vector<Node> N; // Set up nodes
  vector<Arc> E;  // Set up arcs
  vector<int> Head, Tail;

  // lazy constraint
  struct callback_data mydata;
};
