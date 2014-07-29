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

#include "core/encseq_api.h"
#include "core/showtime.h"
#include "core/stack-inlined.h"
#include "core/defined-types.h"
#include "match/esa-seqread.h"
#include "match/esa-smax.h"
#include "match/esa-NE-repeat.h"

typedef struct
{
  GtUword lcp;
  GtUword suf;
  GtUword lb;
  Definedunsignedint left_context;
} Lcp_stackelem;

GT_STACK_DECLARESTRUCT(Lcp_stackelem, 32UL);

int gt_run_NE_repeats_scan(Sequentialsuffixarrayreader *ssar,
                      GtUword searchlength,
                      bool silent,
                      GtProcessNEintervals process_NE_repeat,
                      GT_UNUSED void *process_NE_repeat_data,
                      GtError *err)
{
  bool haserr = false;
  gt_error_check(err);
  GtTimer *nerepeat_progress = NULL;
  GtStackLcp_stackelem lcpstack;
  Lcp_stackelem current_elem;
  GtArrayGtUword suftab_arr;
  GtUword lcp,
          plcp,
          nsuf,
          psuf,
          idx,
          lb,
          pos_corr,
          nonspecials;
  Definedunsignedint bwt_psuf,
                     bwt_nsuf,
                     bwt;

  const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(ssar);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  GT_STACK_INIT(&lcpstack,32UL);
  GT_INITARRAY(&suftab_arr,GtUword);

  gt_showtime_enable();
  if (gt_showtime_enabled())
  {
    nerepeat_progress = gt_timer_new_with_progress_description(
        "finding non-extendible repeats");
  gt_timer_start(nerepeat_progress);
  }

  lcp = 0;
  plcp = 0;
  nsuf = 0;
  psuf = 0;
  lb = 0;
  pos_corr = 0;

  for (idx = 0; idx < nonspecials; idx++)
  {
    if (idx == 0) 
    {
      SSAR_NEXTSEQUENTIALSUFTABVALUE(psuf,ssar);
      current_elem.lcp = plcp;
      current_elem.suf = psuf;
      current_elem.lb = idx;
      if (psuf > 0)
      {
        bwt_psuf.defined = true;
        bwt_psuf = get_left_context(encseq,psuf);
        current_elem.left_context = bwt_psuf;
      } else 
      {
        bwt_psuf.defined = false;
        bwt_psuf.valueunsignedint = 0;
        current_elem.left_context = bwt_psuf;
      }
      GT_STOREINARRAY(&suftab_arr,GtUword,32,psuf);
      GT_STACK_PUSH(&lcpstack,current_elem);
    } else
    {
      SSAR_NEXTSEQUENTIALLCPTABVALUE(lcp,ssar);
      SSAR_NEXTSEQUENTIALSUFTABVALUE(nsuf,ssar);
      GT_STOREINARRAY(&suftab_arr,GtUword,32,nsuf);

      lb = idx-1;

      bwt_nsuf = get_left_context(encseq,nsuf);
      bwt = check_left_context(bwt_psuf, bwt_nsuf);
      bwt_psuf = bwt_nsuf;

      while (GT_STACK_TOP(&lcpstack).lcp > lcp)
      {
        current_elem = GT_STACK_POP(&lcpstack);
        
        if (!current_elem.left_context.defined && current_elem.lcp >= searchlength)
        {
          if (!silent)
          {
            process_NE_repeat(process_NE_repeat_data,
                              encseq,
                              &suftab_arr.spaceGtUword[current_elem.lb-pos_corr],
                              current_elem.lcp,
                              ((idx)-pos_corr-1)-(current_elem.lb-pos_corr-1));
/*            
            printf("LCP: " GT_WU " i: " GT_WU " j: " GT_WU "\n",
                  current_elem.lcp, current_elem.lb, idx-1);
*/
          }
        } 
    
        lb = current_elem.lb;
        GT_STACK_TOP(&lcpstack).left_context = 
                      check_left_context(current_elem.left_context,
                      GT_STACK_TOP(&lcpstack).left_context);
        bwt = check_left_context(current_elem.left_context,bwt);
      }

      if (GT_STACK_TOP(&lcpstack).lcp == lcp)
      {
            GT_STACK_TOP(&lcpstack).left_context =
                                    check_left_context(
                                    GT_STACK_TOP(&lcpstack).left_context,
                                    bwt);
      } else 
      {
        current_elem.lcp = lcp;
        current_elem.suf = psuf;
        current_elem.lb = lb;
        current_elem.left_context = bwt;
        GT_STACK_PUSH(&lcpstack,current_elem);
      }
      if (lcp == 0)
      {
        pos_corr = idx;
//        printf("pos_corr: " GT_WU "\n",pos_corr);
        GT_FREEARRAY(&suftab_arr,GtUword);
        GT_STOREINARRAY(&suftab_arr,GtUword,32,nsuf);
      }
/*
      printf("LCP: " GT_WU "\t",lcp);
      printf("PREVIOUSLCP: " GT_WU "\t",plcp);
      printf("LB: " GT_WU "\t",lb);
      printf("NSUF: " GT_WU "\t",nsuf);
      printf("PSUF: " GT_WU "\n",psuf);
*/
      plcp = lcp;
      psuf = nsuf;
    }
  }
  GT_FREEARRAY(&suftab_arr,GtUword);

  if (nerepeat_progress != NULL)
  {
    gt_timer_show_progress_final(nerepeat_progress, stdout);
    gt_timer_delete(nerepeat_progress);
  }

/*  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
*/
  return haserr ? -1 : 0;
}
