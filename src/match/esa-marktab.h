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

#ifndef ESA_MARKTAB_H
#define ESA_MARKTAB_H

#include "esa-seqread.h"

typedef struct GtESAMarkTab GtESAMarkTab;

GtESAMarkTab* gt_esa_marktab_new(const GtEncseq *encseq);
void gt_esa_marktab_set(GtESAMarkTab *, unsigned int);
bool gt_esa_marktab_get(const GtESAMarkTab *, unsigned int);
void gt_esa_marktab_reset(GtESAMarkTab *);
void gt_esa_marktab_delete(GtESAMarkTab *);

#endif
