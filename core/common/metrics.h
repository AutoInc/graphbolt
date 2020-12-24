#ifndef GRAPHBOLT_METRICS_H_
#define GRAPHBOLT_METRICS_H_
#include <cstdlib>
#include <cstring>
#include <iostream>
namespace pbbs {
inline int parseLine(char* line) {
  // This assumes that a digit will be found and the line ends in " Kb".
  int i = strlen(line);
  const char* p = line;
  while (*p < '0' || *p > '9')
    p++;
  line[i - 3] = '\0';
  i = atoi(p);
  return i;
}

inline double RSSInMB() {  // Note: this value is in KB!
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  if (file == nullptr) {
    return 0;
  }

  while (fgets(line, 128, file) != nullptr) {
    if (strncmp(line, "VmRSS:", 6) == 0) {
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result / 1024.0;
}
}  // namespace pbbs

#endif  // GRAPHBOLT_METRICS_H_
