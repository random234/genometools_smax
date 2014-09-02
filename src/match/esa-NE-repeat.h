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

#ifndef ESA_NE_REPEAT_H
#define ESA_NE_REPEAT_H

#include "esa-seqread.h"
#include "core/encseq_api.h"
#include "core/showtime.h"
#include "match/esa-NE-repeat.h"
/*
typedef void (*GtProcessNEintervals)(void *,
                                     const GtEncseq *,
                                     const GtUword *,
                                     GtUword,
                                     GtUword,
                                     GtUword,
                                     GtUword);
*/
typedef void (*GtProcessNEintervals)(void *,
                                     const GtEncseq *,
                                     const GtUword *,
                                     GtUword,
                                     GtUword);

Definedunsignedint get_lctx(const GtEncseq *,
                                    GtUword);

Definedunsignedint check_lctx(Definedunsignedint,
                                      Definedunsignedint);

bool is_notleftextendible(const GtEncseq *,
                          const GtUword *,
                          const GtUword);
#endif
