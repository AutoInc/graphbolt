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

#include "../core/common/utils.h"
#include "../core/graphBolt/GraphBoltEngine_simple.h"
#include "../core/main.h"
#include <math.h>

#include <utility>

// ======================================================================
// PAGERANKINFO
// ======================================================================
class PageRankInfo {
 public:
  // Should I just use vectors for this?
  uintV n;
  double epsilon;
  double damping;
  long *out_degrees;
  std::shared_ptr<std::map<uintV, double>> answer_base;
  std::shared_ptr<std::map<uintV, double>> answer_inc;

  PageRankInfo()
      : n(0), epsilon(0), damping(0), out_degrees(nullptr), answer_base(
      nullptr), answer_inc(nullptr) {
  }

  PageRankInfo(uintV _n, double _epsilon, double _damping, std::shared_ptr<std::map<uintV, double>> _answer_base,
               std::shared_ptr<std::map<uintV, double>> _answer_inc)
      : n(_n), epsilon(_epsilon), damping(_damping), answer_base(std::move(_answer_base)),
      answer_inc(std::move(_answer_inc)) {
    if (n > 0) {
      out_degrees = newA(long, n);
      parallel_for (uintV i = 0; i < n; i++) { out_degrees[i] = 0; }
    }
  }

  void copy(const PageRankInfo &object) {
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
    epsilon = object.epsilon;
    damping = object.damping;
    answer_base = object.answer_base;
    answer_inc = object.answer_inc;
  }

  ~PageRankInfo() {
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
      writeAdd(&out_degrees[source], (long) 1);
    }
    parallel_for (long i = 0; i < edge_deletions.size; i++) {
      uintV source = edge_deletions.E[i].source;
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
double initial_vertex_value = 0;
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
  v_vertex_value = initial_vertex_value;
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
  vertex_value_next =
      (1 - global_info.damping) / global_info.n  + (global_info.damping * aggregation_value);
}

template<class VertexValueType, class GlobalInfoType>
inline bool isChanged(const VertexValueType &value_curr,
                      const VertexValueType &value_next,
                      GlobalInfoType &global_info) {
  // return (fabs(value_next - value_curr) > global_info.epsilon);
  // comment out this condition because of wrong result.
  return true;
}

// We determine termination based on the diff of sum of vertex values
template<class VertexValueType, class GlobalInfoType>
inline bool isTerminated(const VertexValueType *values_curr,
                         GlobalInfoType &global_info, bool isInc) {
  if (isInc) {
    if (global_info.answer_inc != nullptr) {
      auto&ans = *global_info.answer_inc;
      VertexValueType diff_sum = 0;
      parallel_for (uintV v = 0; v < global_info.n; v++) {
        writeAdd(&diff_sum, fabs(values_curr[v] - ans[v]));
      }
      std::cout << "Inc Diff sum: " << diff_sum;
      return diff_sum < global_info.epsilon;
    }
  } else {
    if (global_info.answer_base != nullptr) {
      auto&ans = *global_info.answer_base;
      VertexValueType diff_sum = 0;
      parallel_for (uintV v = 0; v < global_info.n; v++) {
        writeAdd(&diff_sum, fabs(values_curr[v] - ans[v]));
      }
      std::cout << "Base Diff sum: " << diff_sum;
      return diff_sum < global_info.epsilon;
    }
  }

  return false;
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
      ? (v_value_curr - v_value_prev) / global_info.getOutDegree(v)
      : 0;
}

template<class AggregationValueType, class VertexValueType, class EdgeDataType,
    class GlobalInfoType>
inline bool edgeFunction(const uintV &u, const uintV &v,
                         const EdgeDataType &edge_weight,
                         const VertexValueType &u_value,
                         AggregationValueType &u_change_in_contribution,
                         GlobalInfoType &global_info) {
  return true;
}

// ======================================================================
// INCREMENTAL COMPUTING / DETERMINING FRONTIER
// =====================================================================
template<class GlobalInfoType>
inline void hasSourceChangedByUpdate(const uintV &v, UpdateType update_type,
                                     bool &activateInCurrentIteration,
                                     GlobalInfoType &global_info,
                                     GlobalInfoType &global_info_old) {
  if (global_info.getOutDegree(v) != global_info_old.getOutDegree(v))
    activateInCurrentIteration = true;
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
  double epsilon = config.getOptionDoubleValue("-epsilon", 0.01);
  double damping = config.getOptionDoubleValue("-damping", 0.85);
  int max_iters = config.getOptionLongValue("-maxIters", 10);
  std::string answer_base_path = config.getOptionValue("-answer_base", "");
  std::string answer_inc_path = config.getOptionValue("-answer_inc", "");

  max_iters += 1;

  std::shared_ptr<std::map<uintV, double>> ans_base;
  std::shared_ptr<std::map<uintV, double>> ans_inc;

  if (!answer_base_path.empty()) {
    cout << "Loading answer base file..." << endl;
    ans_base = readAnswer<double>(answer_base_path.c_str());
  }

  if (!answer_inc_path.empty()) {
    cout << "Loading answer inc file..." << endl;
    ans_inc = readAnswer<double>(answer_inc_path.c_str());
  }

  PageRankInfo
      global_info
      (n, epsilon, damping, ans_base,
       ans_inc);
  parallel_for (uintV i = 0; i < n; i++) {
    global_info.out_degrees[i] = G.V[i].getOutDegree();
  }

  cout << "Initializing engine ....\n";
  GraphBoltEngineSimple<vertex, double, double, PageRankInfo> engine(
      G, max_iters, global_info, false, config);
  engine.init();
  cout << "Finished initializing engine\n";

  engine.run();
}
