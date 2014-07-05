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
#include "match/esa-seqread.h"
#include "match/esa-smax.h"
#include "core/defined-types.h"

/* vorraussichtlice Stackgroesse muss geschaetzt werden um
 * konstante Stackgroesse sinnvoll zu waehlen.
 */

typedef struct
{
  GtUword lcp;
  GtUword lb;
  GtUword suf;
  //Definedunsignedint left_context;
  //char left_context;
  //bool left_extendible;
} Lcp_stackelem;

GT_STACK_DECLARESTRUCT(Lcp_stackelem, 32UL);

int gt_run_NE_repeats(Sequentialsuffixarrayreader *ssar,
                      GT_UNUSED GtUword searchlength,
                      GT_UNUSED bool silent,
                      GT_UNUSED GtProcessSmaxpairs process_smaxpairs,
                      GT_UNUSED void *process_smaxpairsdata,
                      GtError *err)
{
  bool haserr = false;
  gt_error_check(err);
  GtTimer *nerepeat_progress = NULL;
  GtStackLcp_stackelem lcpstack;
  Lcp_stackelem current_elem;
//  const ESASuffixptr *suftab;
  bool interval;
  GT_UNUSED GtUword lcp,
          plcp,
          psuf,
          ppsuf,
          lb,
          plb,
          idx,
          nonspecials,
          seqnum,
          clmax;

  const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(ssar);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
//  suftab = gt_suffixarraySequentialsuffixarrayreader(ssar)->suftab;
  GT_STACK_INIT(&lcpstack,32UL);

/*  gt_showtime_enable();
  if (gt_showtime_enabled())
  {
    nerepeat_progress = gt_timer_new_with_progress_description(
        "finding non-extendible repeats");
  gt_timer_start(nerepeat_progress);
  }
*/

  clmax = 0;
  seqnum = 0;
  interval = false;
  lcp = 0;
  plcp = 0;
  lb = 0;
  plb = 0;
  psuf = 0;
  ppsuf = 0;

  current_elem.lcp = lcp;
  current_elem.lb = lb;
  current_elem.suf = psuf;

  GT_STACK_PUSH(&lcpstack,current_elem);

  for (idx = 0; idx < nonspecials; idx++)
  {
    SSAR_NEXTSEQUENTIALLCPTABVALUE(lcp,ssar);
    SSAR_NEXTSEQUENTIALSUFTABVALUE(psuf,ssar);
    lb = idx;
    printf("\nLCP: " GT_WU "\t",lcp);
    printf("LB: " GT_WU "\t",lb);
    printf("SUF: " GT_WU "\n",psuf);

    if (interval)
    {
      bool left_context = false;
      printf("in interval\n");
      while (GT_STACK_TOP(&lcpstack).lcp > plcp)
      {
        char left_char1 = 0;
        char left_char2 = 0;
        current_elem = GT_STACK_POP(&lcpstack);
        /*
        printf("" GT_WU " " GT_WU " " GT_WU "\n",current_elem.lcp,
                                                current_elem.lb,
                                                current_elem.suf);
        */
        if (left_context)
        {
        } else
        {
          if (current_elem.suf > 0) 
          {
            left_char1 = gt_encseq_get_decoded_char(encseq,
                                                    current_elem.suf-1,
                                                    GT_READMODE_FORWARD);
          } else
          {
            left_char1 = '$';
          }
          if (ppsuf > 0)  
          { 
            left_char2 = gt_encseq_get_decoded_char(encseq,
                                                    ppsuf-1,
                                                    GT_READMODE_FORWARD);
          } else
          {
            left_char2 = '$';
          }

          if (ISSPECIAL(left_char1) || left_char1 != left_char2)
          {
            printf("" GT_WU " " GT_WU " " GT_WU " F ",current_elem.lcp,
                                                    seqnum,
                                                    current_elem.suf);
            printf("" GT_WU " " GT_WU " " GT_WU "\n",current_elem.lcp,
                                                    seqnum,
                                                    ppsuf);
            left_context = true;
          }
          left_char1 = 0;
          left_char2 = 0;
                    
        }
        // SAVE LAST POP
      }
      interval = false;
      /*
      printf("" GT_WU " ",plcp);
      printf("" GT_WU " ",plb);
      printf("" GT_WU "\n",ppsuf);
      */
      clmax = plcp;
    }

    if (clmax == lcp)
    { 
      printf("clmax == lcp\n");
      current_elem.lcp = lcp;
      current_elem.lb = lb;
      current_elem.suf = psuf;
      GT_STACK_PUSH(&lcpstack,current_elem);
    }

    if (clmax > lcp)
    {
      printf("clmax > lcp\n");
      interval = true;
      plcp = lcp;
      plb = lb;
      ppsuf = psuf;
      clmax = lcp;
    }

    if (clmax < lcp)
    {
      printf("clmax < lcp \t interval=true\n");
      current_elem.lcp = lcp;
      current_elem.lb = lb;
      current_elem.suf = psuf;
      GT_STACK_PUSH(&lcpstack,current_elem);
      clmax = lcp;
    }

/*    
      GT_STACK_PUSH(&lcpstack,current_elem);
      SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
      SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
      GT_STACK_TOP(&lcpstack).lcp
      current_elem = GT_STACK_POP(&lcpstack);
*/
  }

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
