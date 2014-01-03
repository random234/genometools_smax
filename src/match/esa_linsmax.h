#ifndef ESA_LINSMAX_H
#define ESA_LINSMAX_H

#include "core/logger_api.h"
#include "core/error_api.h"
#include "match/esa-seqread.h"

inline int gt_runlinsmax(GtStrArray *,GtUword, bool, bool, bool, GtLogger *, GtError *);
inline bool gt_linsmax_verify_supmax(Sequentialsuffixarrayreader *, GtUword, GtUword);

#endif
