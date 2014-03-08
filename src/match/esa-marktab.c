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
