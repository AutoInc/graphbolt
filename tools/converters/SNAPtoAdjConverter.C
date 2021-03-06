// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
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

// Converts a SNAP graph (http://snap.stanford.edu/data/index.html) to
// Ligra adjacency graph format. To symmetrize the graph, pass the "-s"
// flag. For undirected graphs on SNAP, the "-s" flag must be passed
// since each edge appears in only one direction

#include "../../core/common/parallel.h"
#include "../../core/common/parseCommandLine.h"
#include "../common/graphIO.h"

int parallel_main(int argc, char *argv[]) {
  commandLine P(argc, argv, "[-s] <input SNAP file> <output Ligra file>");
  char *iFile = P.getArgument(1);
  char *oFile = P.getArgument(0);
  bool sym = P.getOption("-s");
  bool weighted = P.getOption("-w");
  cout << "Reading graph and creating \n";
  edgeArray G;
  if (weighted) {
    wghEdgeArray G = readWghSNAP(iFile);
    cout << "Writing to output file\n";
    writeWghGraphToFile(wghGraphFromWghEdges(G, sym), oFile);
  } else {
    edgeArray G = readSNAP(iFile);
    cout << "Writing to output file\n";
    writeGraphToFile(graphFromEdges(G, sym), oFile);
  }
}
