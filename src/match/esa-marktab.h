#ifndef ESA_MARKTAB_H
#define ESA_MARKTAB_H

#include "esa-seqread.h"

typedef struct GtESAMarkTab GtESAMarkTab;

GtESAMarkTab* gt_esa_marktab_new(const GtEncseq *encseq);
void gt_esa_marktab_set(GtESAMarkTab *, int, bool);
int gt_esa_marktab_get(GtESAMarkTab *, int);
int gt_esa_marktab_size(GtESAMarkTab *tab);
void gt_esa_marktab_reset(GtESAMarkTab *);

#endif


