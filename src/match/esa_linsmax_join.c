#include "esa-seqread.h"
#include "core/encseq_api.h"
#include "core/str_api.h"
#include "core/showtime.h"
#include "match/esa_linsmax_join.h"


GtStr* gt_linsmax_join(GtStrArray *index_arr, GT_UNUSED GtLogger *logger, GtError *err)
{
	gt_assert(index_arr);
  GtStr *indexname = gt_str_new();
	GtAlphabet *alpha = gt_alphabet_new_from_file(gt_str_array_get(index_arr,0), err);
	GtEncseqBuilder *eb = gt_encseq_builder_new(alpha);
	GT_UNUSED GtEncseqReader *er;
	GT_UNUSED GtEncseq *encseq;
	GT_UNUSED Sequentialsuffixarrayreader *sarr;
	GT_UNUSED GtUword i,j = 0;

	gt_encseq_builder_enable_multiseq_support(eb);
	gt_encseq_builder_create_esq_tab(eb);
	gt_encseq_builder_create_des_tab(eb);
	gt_encseq_builder_create_ssp_tab(eb);
	gt_encseq_builder_create_sds_tab(eb);



	for(i=0;i<=gt_str_array_size(index_arr);i++)
	{
		


/*		sarr = gt_newSequentialsuffixarrayreaderfromfile(gt_str_array_get(index_arr,i),			                                                   SARR_LCPTAB |
		 SARR_SUFTAB | 	                                                    SARR_ESQTAB,									                                     SEQ_mappedboth,								                                    //SEQ_scan,
     logger,
     err);
		for(i=0;i<=gt_
		gt_encseq_builder_add_encoded(eb, const GtUchar *str, GtUword strlen, const char *desc);
	*/
	}
	return indexname;
}
