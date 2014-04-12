/*
  Copyright (c) 2014 Philipp Lutz Carpus <random234@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef ESA_SMAX_MAP_H
#define ESA_SMAX_MAP_H

#include "core/str_api.h"
#include "match/esa_visitor.h"
#include "match/esa-seqread.h"
#include "match/esa-smax.h"

typedef struct GtESASmaxLcpintervalsVisitor GtESASmaxLcpintervalsVisitor;

const GtESAVisitorClass* gt_esa_smax_lcpitvs_visitor_class(void);
GtESAVisitor* gt_esa_smax_lcpitvs_visitor_new(Sequentialsuffixarrayreader *,
                                              GtUword,
                                              bool,
                                              bool,
                                              GtProcessSmaxpairs,
                                              void *,
                                              GtTimer *);
bool gt_esa_smax_lcpitvs_visitor_get_info(GtESAVisitorInfo *);
void gt_esa_smax_lcpitvs_visitor_set_info(GtESAVisitorInfo *, bool);

#endif
