#include "esa-seqread.h"
#include "esa-marktab.h"

struct GtESAMarkTab {
	int size;
	bool *alpha;
};

GtESAMarkTab* gt_esa_marktab_new(const GtEncseq *encseq)
{
	GtESAMarkTab *tab;
	tab = gt_malloc(sizeof (GtESAMarkTab));
	tab->size = gt_alphabet_size(gt_encseq_alphabet(encseq));
	tab->alpha = gt_malloc(sizeof (bool) * tab->size);
	return tab;
}

void gt_esa_marktab_set(GtESAMarkTab *tab, int pos, bool val)
{
	tab->alpha[pos] = val;
}

int gt_esa_marktab_get(GtESAMarkTab *tab, int pos)
{
	return tab->alpha[pos];
}

int gt_esa_marktab_size(GtESAMarkTab *tab)
{
	return tab->size;
}

void gt_esa_marktab_reset(GtESAMarkTab *tab)
{
	int i;
	for(i=0;i<=tab->size;i++)
	{
		tab->alpha[i] =false;
	}
}
void gt_esa_marktab_delete(GtESAMarkTab *tab)
{
	gt_free(tab);
}
