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

#include "esa-seqread.h"
#include "esa-marktab.h"

struct GtESAMarkTab {
  unsigned int size;
  bool *alpha;
};

GtESAMarkTab* gt_esa_marktab_new(const GtEncseq *encseq)
{
  GtESAMarkTab *tab = gt_malloc(sizeof (*tab));
  tab->size = gt_alphabet_size(gt_encseq_alphabet(encseq));
  tab->alpha = gt_malloc(sizeof (*tab->alpha) * tab->size);
  return tab;
}

void gt_esa_marktab_set(GtESAMarkTab *tab, unsigned int pos)
{
  tab->alpha[pos] = true;
}

bool gt_esa_marktab_get(const GtESAMarkTab *tab, unsigned int pos)
{
  return tab->alpha[pos];
}

void gt_esa_marktab_reset(GtESAMarkTab *tab)
{
  unsigned int idx;

  for (idx=0;idx<tab->size;idx++)
  {
    tab->alpha[idx] =false;
  }
}

void gt_esa_marktab_delete(GtESAMarkTab *tab)
{
  if (tab != NULL)
  {
    gt_free(tab->alpha);
    gt_free(tab);
  }
}
