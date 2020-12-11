// Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// NOTE: The edge data type header file should then be included as the first header
// file at the top of the user program.
#include "PHP_edgeData.h"

#include "../core/common/utils.h"
#include "../core/graphBolt/GraphBoltEngine_simple.h"
//#include "../core/graphBolt/GraphBoltEngine_complex.h"
#include "../core/main.h"
#include <math.h>

// ======================================================================
// PHPINFO
// ======================================================================
class PHPInfo {
 public:
  uintV n;
  uintV source;
  double epsilon;
  double damping;
  long *out_degrees;

  PHPInfo() : n(0), source(0), epsilon(0), damping(0), out_degrees(nullptr) {
  }

  PHPInfo(uintV _n, uintV source, double _epsilon, double _damping)
      : n(_n), source(source), epsilon(_epsilon), damping(_damping) {
    if (n > 0) {
      out_degrees = newA(long, n);
      parallel_for (uintV i = 0; i < n; i++) { out_degrees[i] = 0; }
    }
  }

  void copy(const PHPInfo &object) {
    if (object.n > n) {
      if (n == 0) {
        n = object.n;
        out_degrees = newA(long, n);
      } else {
        // realloc
        n = object.n;
        out_degrees = renewA(long, out_degrees, n);
      }
    }
    long min_n = std::min(object.n, n);
    parallel_for (uintV i = 0; i < min_n; i++) {
      out_degrees[i] = object.out_degrees[i];
    }
    source = object.source;
    epsilon = object.epsilon;
    damping = object.damping;
  }

  ~PHPInfo() {
    if (n > 0)
      deleteA(out_degrees);
  }

  void processUpdates(edgeArray &edge_additions, edgeArray &edge_deletions) {
    // Increase out_degrees array size
    if (edge_additions.maxVertex >= n) {
      uintV n_old = n;
      n = edge_additions.maxVertex + 1;
      out_degrees = renewA(long, out_degrees, n);
      parallel_for (uintV i = n_old; i < n; i++) { out_degrees[i] = 0; }
    }

    parallel_for (long i = 0; i < edge_additions.size; i++) {
      uintV source = edge_additions.E[i].source;
      uintV destination = edge_additions.E[i].destination;
      writeAdd(&out_degrees[source], (long) 1);
    }
    parallel_for (long i = 0; i < edge_deletions.size; i++) {
      uintV source = edge_deletions.E[i].source;
      uintV destination = edge_deletions.E[i].destination;
      writeAdd(&out_degrees[source], (long) -1);
    }
  }

  void cleanup() {
  }

  inline long getOutDegree(const uintV &v) {
    return (v < n) ? out_degrees[v] : 0;
  }
};

// ======================================================================
// AGGREGATEVALUE AND VERTEXVALUE INITIALIZATION
// ======================================================================
double initial_aggregation_value = 0;
double aggregation_value_identity = 0;
double vertex_value_identity = 0;
template<class AggregationValueType, class GlobalInfoType>
inline void
initializeAggregationValue(const uintV &v,
                           AggregationValueType &v_aggregation_value,
                           const GlobalInfoType &global_info) {
  v_aggregation_value = initial_aggregation_value;
}

template<class VertexValueType, class GlobalInfoType>
inline void initializeVertexValue(const uintV &v,
                                  VertexValueType &v_vertex_value,
                                  const GlobalInfoType &global_info) {
  v_vertex_value = (v == global_info.source ? 1 : 0);
}

template<class AggregationValueType>
inline AggregationValueType &aggregationValueIdentity() {
  return aggregation_value_identity;
}

template<class VertexValueType>
inline VertexValueType &vertexValueIdentity() {
  return vertex_value_identity;
}
// ======================================================================
// ACTIVATE VERTEX/COMPUTE VERTEX FOR A GIVEN ITERATION
// ======================================================================
template<class GlobalInfoType>
inline bool forceActivateVertexForIteration(const uintV &v, int iter,
                                            const GlobalInfoType &global_info) {
  if (iter == 1) {
    return true;
  } else {
    return false;
  }
}

template<class GlobalInfoType>
inline bool forceComputeVertexForIteration(const uintV &v, int iter,
                                           const GlobalInfoType &global_info) {
  if (iter == 1) {
    return true;
  } else {
    return false;
  }
}
inline bool shouldUseDelta(int iter) {
  if (iter == 1) {
    return false;
  } else {
    return true;
  }
}

// ======================================================================
// ADD TO OR REMOVE FROM AGGREGATION VALUES
// ======================================================================
template<class AggregationValueType, class GlobalInfoType>
inline void addToAggregation(const AggregationValueType &incoming_value,
                             AggregationValueType &aggregate_value,
                             GlobalInfoType &global_info) {
  aggregate_value += incoming_value;
}

template<class AggregationValueType, class GlobalInfoType>
inline void addToAggregationAtomic(const AggregationValueType &incoming_value,
                                   AggregationValueType &aggregate_value,
                                   GlobalInfoType &global_info) {
  writeAdd(&aggregate_value, incoming_value);
}

template<class AggregationValueType, class GlobalInfoType>
inline void removeFromAggregation(const AggregationValueType &incoming_value,
                                  AggregationValueType &aggregate_value,
                                  GlobalInfoType &global_info) {
  aggregate_value -= incoming_value;
}

template<class AggregationValueType, class GlobalInfoType>
inline void
removeFromAggregationAtomic(const AggregationValueType &incoming_value,
                            AggregationValueType &aggregate_value,
                            GlobalInfoType &global_info) {
  writeAdd(&aggregate_value, -incoming_value);
}

// ======================================================================
// VERTEX COMPUTE FUNCTION AND DETERMINE END OF COMPUTATION
// ======================================================================
template<class AggregationValueType, class VertexValueType,
    class GlobalInfoType>
inline void computeFunction(const uintV &v,
                            const AggregationValueType &aggregation_value,
                            const VertexValueType &vertex_value_curr,
                            VertexValueType &vertex_value_next,
                            GlobalInfoType &global_info) {
  if (v == global_info.source == v) {
    vertex_value_next = 1;
  } else {
    vertex_value_next =
        global_info.damping * aggregation_value;
  }
}

template<class VertexValueType, class GlobalInfoType>
inline bool isChanged(const VertexValueType &value_curr,
                      const VertexValueType &value_next,
                      GlobalInfoType &global_info) {
  // return (fabs(value_next - value_curr) > global_info.epsilon);
  return true;
}

// ======================================================================
// EDGE FUNCTIONS
// ======================================================================
template<class AggregationValueType, class VertexValueType,
    class StaticInfoType>
inline void sourceChangeInContribution(
    const uintV &v, AggregationValueType &v_change_in_contribution,
    const VertexValueType &v_value_prev, const VertexValueType &v_value_curr,
    StaticInfoType &global_info) {
  v_change_in_contribution =
      (global_info.getOutDegree(v) != 0)
      ? v_value_curr - v_value_prev
      : 0;
}

template<class AggregationValueType, class VertexValueType, class EdgeDataType,
    class GlobalInfoType>
inline bool edgeFunction(const uintV &u, const uintV &v,
                         const EdgeDataType &edge_weight,
                         const VertexValueType &u_value,
                         AggregationValueType &u_change_in_contribution,
                         GlobalInfoType &global_info) {

//  if(v == global_info.source) {
//    return false;
//  } else {
    u_change_in_contribution *= edge_weight.weight;
    return true;
//  }
}

//template<class AggregationValueType, class VertexValueType, class EdgeDataType,
//    class GlobalInfoType>
//inline bool edgeFunctionDelta(const uintV &u, const uintV &v,
//                              const EdgeDataType &edge_weight,
//                              const VertexValueType &u_value_prev,
//                              const VertexValueType &u_value_curr,
//                              AggregationValueType &u_change_in_contribution,
//                              GlobalInfoType &global_info) {
//  if (v == global_info.source) {
//    return false;
//  } else {
//    u_change_in_contribution =
//        (u_value_curr - u_value_prev) * edge_weight.weight;
//    return true;
//  }
//}

// ======================================================================
// INCREMENTAL COMPUTING / DETERMINING FRONTIER
// =====================================================================
template<class GlobalInfoType>
inline void hasSourceChangedByUpdate(const uintV &v, UpdateType update_type,
                                     bool &activateInCurrentIteration,
                                     GlobalInfoType &global_info,
                                     GlobalInfoType &global_info_old) {
  // TODO: need this?
//  if (global_info.getOutDegree(v) != global_info_old.getOutDegree(v))
//    activateInCurrentIteration = true;
}

template<class GlobalInfoType>
inline void hasDestinationChangedByUpdate(const uintV &v,
                                          UpdateType update_type,
                                          bool &activateInCurrentIteration,
                                          GlobalInfoType &global_info,
                                          GlobalInfoType &global_info_old) {
}

// ======================================================================
// HELPER FUNCTIONS
// ======================================================================
template<class AggregationValueType, class VertexValueType,
    class GlobalInfoType>
void printHistory(const uintV &v, AggregationValueType **agg_values,
                  VertexValueType **vertex_values, GlobalInfoType &info,
                  int history_iterations) {
  for (int iter = 0; iter < history_iterations; iter++) {
    cout << iter << "," << agg_values[iter][v] << "," << vertex_values[iter][v]
         << "\n";
  }
}

template<class GlobalInfoType>
void printAdditionalData(ofstream &output_file, const uintV &v,
                         GlobalInfoType &info) {
}

// ======================================================================
// COMPUTE FUNCTION
// ======================================================================
template<class vertex>
void compute(graph<vertex> &G, commandLine config) {
  uintV n = G.n;
  int max_iters = config.getOptionLongValue("-maxIters", 10);
  uintV source = config.getOptionLongValue("-source", 0);
  max_iters += 1;
  double epsilon = 0.01;
  double damping = 0.8;

  PHPInfo global_info(n, source, epsilon, damping);
  parallel_for (uintV i = 0; i < n; i++) {
    global_info.out_degrees[i] = G.V[i].getOutDegree();
  }

  cout << "Initializing engine ....\n";
  GraphBoltEngineSimple<vertex, double, double, PHPInfo> engine(
      G, max_iters, global_info, false, config);
  engine.init();
  cout << "Finished initializing engine\n";

  engine.run();
}
