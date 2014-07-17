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

#ifndef ESA_SMAX_H
#define ESA_SMAX_H

#include "core/encseq_api.h"
#include "core/types_api.h"
#include "core/arraydef.h"
#include "core/bittab_api.h"

typedef void (*GtProcessSmaxpairs)(void *,
                                  const struct GtEncseq *,
                                  const GtUword,
                                  const GtUword *,
                                  const GtUword);

typedef void (*GtProcessNEintervals)(void *,
                                    const struct GtEncseq *,
                                    const GtUword *,
                                    GtUword,
                                    GtUword,
                                    GtUword,
                                    GtUword);

bool gt_esa_smax_verify_supmax(const GtEncseq *,
                              const GtUword *,
                              const GtUword,
                              GtBittab *);
bool gt_esa_smax_verify_supmax_count(const GtEncseq *,
                              const GtUword *,
                              const GtUword,
                              GtBittab *,
                              GtUword *);
#endif
