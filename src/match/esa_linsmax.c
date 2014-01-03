#include "esa-seqread.h"
#include "core/encseq_api.h"
#include "core/str_api.h"
#include "core/showtime.h"
#include "match/esa_linsmax.h"


int gt_runlinsmax(GtStrArray *inputindex,
            GT_UNUSED GtUword searchlength,
            GT_UNUSED bool absolute,
            GT_UNUSED bool silent,
            GT_UNUSED bool outedges,
            GtLogger *logger,
            GtError *err)
{
gt_assert(inputindex);
bool haserr = false;
GT_UNUSED Sequentialsuffixarrayreader *ssar;
GtTimer *linsmaxprogress = NULL;
gt_error_check(err);
ssar = gt_newSequentialsuffixarrayreaderfromfile(gt_str_array_get(inputindex,0),
                                                 SARR_LCPTAB |
                                                 SARR_SUFTAB |
                                                 SARR_ESQTAB,
                                                 false,
                                                 //SEQ_scan,
                                                 logger,
                                                 err);
const GtUword *suftab;
const Suffixarray *sa;
GT_UNUSED char method = 'D';
GT_UNUSED int pos_corr_t,pos_corr_s = 0;
GtUword lcpvalue, 
							previoussuffix = 0,
							idx,
							nonspecials,
							sumsuftab = 0,
							sumlcptab = 0,
							max = 0,
							lb = 0,
							rb = 0,
							s,
							t;

sa = gt_suffixarraySequentialsuffixarrayreader(ssar);
suftab = gt_suftabSequentialsuffixarrayreader(ssar);

gt_showtime_enable();
if (gt_showtime_enabled())
{
	linsmaxprogress = gt_timer_new_with_progress_description("attempting to find supermaximal repeats");
	gt_timer_start(linsmaxprogress);
}
nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);

//printf("nonspecials: %lu \n",nonspecials);

for (idx = 0; idx < nonspecials; idx++)
{
	
	SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
	sumlcptab += lcpvalue; /* silly but guarantees that loop is not eliminated by compiler */
	SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
//	printf("SUF: %lu LCP: %lu\n",previoussuffix,lcpvalue);
	sumsuftab += previoussuffix; /* silly but guarantees that loop is not eliminated by compiler */

	if(lcpvalue > max) 
	{
		max = lcpvalue;
		rb = idx+1;
	} else
	{
//		printf("End of local maxima lb: %lu rb: %lu\n",lb,rb);

//gt_linsmax_verify_supmax(ssar,lb,rb)

  if(gt_linsmax_verify_supmax(ssar,lb,rb))
	{		
		for(s=lb;s<rb;s++)
		{
			for(t=s+1;t<=rb;t++)
			{
				if(suftab[s] < suftab[t])
        {																																			             
          if(absolute)
          {
            if(!silent)
            {
              printf("%2lu %3lu %3c %2lu %2lu\n",max, suftab[s], method, max, suftab[t]);
            }
          }
          else
          {
            pos_corr_t = gt_encseq_seqstartpos(sa->encseq, gt_encseq_seqnum(sa->encseq,suftab[t]));
            pos_corr_s = gt_encseq_seqstartpos(sa->encseq, gt_encseq_seqnum(sa->encseq,suftab[s]));
            if(!silent)
            {
              printf("%2lu %4lu %3lu %3c %2lu %4lu %2lu\n",max, gt_encseq_seqnum(sa->encseq,suftab[s]), suftab[s]-pos_corr_s, method, max, gt_encseq_seqnum(sa->encseq,suftab[t]), suftab[t]-pos_corr_t);
            }
          }
        } else 
				{
					if(absolute)
          {
	          if(!silent)
	          {
              printf("%2lu %3lu %3c %2lu %2lu\n",max, suftab[t], method, max, suftab[s]);
            }
          }
          else
	        {			
            pos_corr_t = gt_encseq_seqstartpos(sa->encseq, gt_encseq_seqnum(sa->encseq,suftab[t]));
            pos_corr_s = gt_encseq_seqstartpos(sa->encseq, gt_encseq_seqnum(sa->encseq,suftab[s]));
            if(!silent)
            {
              printf("%2lu %4lu %3lu %3c %2lu %4lu %2lu\n",max, gt_encseq_seqnum(sa->encseq,suftab[t]), suftab[t]-pos_corr_t, method, max, gt_encseq_seqnum(sa->encseq,suftab[s]), suftab[s]-pos_corr_s);
	          }
          }
        }
      }
			lb = idx+1;
			max=0;
		}
	}
}
}

if (linsmaxprogress != NULL)
{
	//gt_timer_show_progress(smaxprogress,"%ld.%06lds real %lds user %lds system", stdout);
	gt_timer_show_progress_final(linsmaxprogress, stdout);
	gt_timer_delete(linsmaxprogress);
}

if (ssar != NULL)
{
	gt_freeSequentialsuffixarrayreader(&ssar);
}

return haserr ? -1 : 0;
}

bool gt_linsmax_verify_supmax(Sequentialsuffixarrayreader *ssar, GtUword lb, GtUword rb) 
{
	const GtUword *suftab;
	const Suffixarray *sa;
	bool marktab[GT_DNAALPHASIZE+1];
	int i = 0;

	sa = gt_suffixarraySequentialsuffixarrayreader(ssar);
	suftab = gt_suftabSequentialsuffixarrayreader(ssar);

	for(i=0;i<=GT_DNAALPHASIZE;i++)
  {
	  marktab[i] = false;
	}
//	printf("lb: %lu rb: %lu\n",lb,rb);
	for(i=lb;i<=rb;i++)
	{
//		printf("%d",i);
//		printf("%lu",suftab[i]);
		if(!suftab[i]==0)
		{
			if(ISNOTSPECIAL(gt_encseq_get_encoded_char(sa->encseq,suftab[i]-1,0)))
			{
				
				if(marktab[gt_encseq_get_encoded_char(sa->encseq,suftab[i]-1,0)] == true)
				{
					return false;
        }
        marktab[gt_encseq_get_encoded_char(sa->encseq,suftab[i]-1,0)] = true;
      }
    }
  }
  return true;
}
