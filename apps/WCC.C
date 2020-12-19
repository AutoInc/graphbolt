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

#ifdef EDGEDATA
// NOTE: The edge data type header file should then be included as the first header
// file at the top of the user program.
// #include "WCC_edgeData.h"
#endif

#include "../core/common/utils.h"
#include "../core/graphBolt/KickStarterEngine.h"
#include "../core/main.h"
#include <math.h>

#define MAX_DISTANCE 65535

// ======================================================================
// SSSPINFO
// ======================================================================
class WccInfo {
public:

  WccInfo(){}

  void copy(const WccInfo &object) {}

  void processUpdates(edgeArray &edge_additions, edgeArray &edge_deletions) {}

  void cleanup() {}
};

// ======================================================================
// VERTEXVALUE INITIALIZATION
// ======================================================================
template <class VertexValueType, class GlobalInfoType>
inline void initializeVertexValue(const uintV &v,
                                  VertexValueType &v_vertex_value,
                                  const GlobalInfoType &global_info) {
  v_vertex_value = v;
}

// ======================================================================
// ACTIVATE VERTEX/COMPUTE VERTEX FOR FIRST ITERATION
// ======================================================================
template <class GlobalInfoType>
inline bool frontierVertex(const uintV &v, const GlobalInfoType &global_info) {
    return true;
}

// ======================================================================
// EDGE FUNCTION
// ======================================================================
template <class VertexValueType, class EdgeDataType, class GlobalInfoType>
inline bool
edgeFunction(const uintV &u, const uintV &v, const EdgeDataType &edge_weight,
             const VertexValueType &u_value, VertexValueType &v_value,
             GlobalInfoType &global_info) {
  cout << u << ":" << u_value << "->" << v << ":" << v_value << endl;
  if (u_value == MAX_DISTANCE) {
    return false;
  } else {
    VertexValueType old_v = v_value;
#ifdef EDGEDATA
    // v_value = min(u_value, v_value);
    v_value = u_value;
#else
    // v_value = min(u_value, v_value);
    v_value = u_value;
#endif
    return true;
  }
}

// ======================================================================
// SHOULDPROPAGATE
// ======================================================================
// shouldPropagate condition for deciding if the value change in
// updated graph violates monotonicity
template <class VertexValueType, class GlobalInfoType>
inline bool shouldPropagate(const VertexValueType &old_value,
                            const VertexValueType &new_value,
                            GlobalInfoType &global_info) {
  return (new_value > old_value);
}

// ======================================================================
// HELPER FUNCTIONS
// ======================================================================
template <class GlobalInfoType>
void printAdditionalData(ofstream &output_file, const uintV &v,
                         GlobalInfoType &info) {}

// ======================================================================
// COMPUTE FUNCTION
// ======================================================================
template <class vertex> void compute(graph<vertex> &G, commandLine config) {
  WccInfo global_info;

  cout << "Initializing engine ....\n";
  KickStarterEngine<vertex, uint16_t, WccInfo> engine(G, global_info, config);
  engine.init();
  cout << "Finished initializing engine\n";
  engine.run();
}
