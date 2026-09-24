#ifndef XCORR_POSTPROCESS_H
#define XCORR_POSTPROCESS_H
#include <vector>
struct st_xcorr;
class OutputTransaction;
void run_postprocess(const st_xcorr &xc, const char *master_prm,
                     const std::vector<int> &xpos, const std::vector<int> &ypos,
                     OutputTransaction &outputs);
#endif
