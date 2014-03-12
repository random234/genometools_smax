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


