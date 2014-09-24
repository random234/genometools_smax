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
  GtUword lb;
  GtUword lcp;
} Lcp_stackelem;

GT_STACK_DECLARESTRUCT(Lcp_stackelem, 32UL);

bool check_NLE(const GtEncseq *encseq, const GtUword *suftab, GtUword suftab_size)
{
  GtUword s;
  GtUword suf;
  GtUchar cc,
          cc_comp;
  suf = suftab[0];
//  printf("SUF: " GT_WU "\n",suf);
  if (suf != 0)
  {
    cc_comp = gt_encseq_get_encoded_char(encseq,suf-1,
                                         GT_READMODE_FORWARD);
    if (ISSPECIAL(cc_comp))
    {
      return true;
    }
  } else 
  {
    return true;
  }
//  printf("cc_comp: %d \n", cc_comp);
//  printf("suftab_size: " GT_WU "\n",suftab_size);
  for (s = 1; s <= suftab_size;s++)
  {
    suf = suftab[s];
//    printf("SUFcc: " GT_WU "\n",suf);
    cc = gt_encseq_get_encoded_char(encseq,suf-1,
                                    GT_READMODE_FORWARD);
  //  printf("cc: %d \n", cc_comp);
    if (ISSPECIAL(cc))
    {
      return true;
    } else 
    {
      if (cc_comp != cc)
      {
        return true;
      }
    }
  }
  return false;
}

int gt_run_NE_repeats_bioscan(Sequentialsuffixarrayreader *ssar,
                      GT_UNUSED GtUword searchlength,
                      bool silent,
                      GT_UNUSED GtProcessNEintervals process_NE_repeat,
                      GT_UNUSED void *process_NE_repeat_data,
                      GtError *err)
{
  bool haserr = false;
  gt_error_check(err);
  GtTimer *nerepeat_progress = NULL;
  GtStackLcp_stackelem lcpstack;
  Lcp_stackelem currelem;
  GtArrayGtUword suftab_arr;
  GtUword nonspecials,
          lcp,
          plcp,
          psuf,
          prevNE,
          lb,
          idx,
          poscorr_flush = 0,
          poscorr_mem = 0;

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
 
  lb = 0;
  prevNE = 0;
  plcp = -1;
  
  currelem.lb = 0;
  currelem.lcp = 0;
  GT_STACK_PUSH(&lcpstack,currelem);
  for (idx = 0; idx < nonspecials; idx++)
  {
    SSAR_NEXTSEQUENTIALLCPTABVALUE(lcp,ssar);
    SSAR_NEXTSEQUENTIALSUFTABVALUE(psuf,ssar);
    printf("SUF: " GT_WU " LCP: " GT_WU "\n",psuf,lcp);

    if (lcp >= searchlength)
    {
      if (lcp == plcp)
      {
        GT_STOREINARRAY(&suftab_arr,GtUword,32,psuf);
      }
      
      if (lcp > plcp)
      {
        currelem.lb = idx;
        currelem.lcp = lcp;
        GT_STACK_PUSH(&lcpstack,currelem);
        GT_STOREINARRAY(&suftab_arr,GtUword,32,psuf);
      } 
    } 

    if (lcp < plcp && plcp >= searchlength)
    {
      GT_STOREINARRAY(&suftab_arr,GtUword,32,psuf);
      while (!(GT_STACK_TOP(&lcpstack).lcp <= lcp))
      {
        if (GT_STACK_TOP(&lcpstack).lcp > 0)
        {
          currelem = GT_STACK_POP(&lcpstack);
          lb = currelem.lb;
          if (prevNE >= currelem.lb)
          { 
            if (!silent)
            {
              process_NE_repeat(process_NE_repeat_data,
                              encseq,
                              &suftab_arr.spaceGtUword[currelem.lb-
                              (poscorr_flush+poscorr_mem)],
                              currelem.lcp,
                              ((idx+1)-(poscorr_flush+poscorr_mem))-
                              (currelem.lb-(poscorr_flush+poscorr_mem)));
            }
          } else
          {
            if (check_NLE(encseq, &suftab_arr.spaceGtUword[currelem.lb-
                        poscorr_flush+poscorr_mem],
                        (idx-(poscorr_flush+poscorr_mem))-
                        (currelem.lb-(poscorr_flush+poscorr_mem))))
            {
              prevNE = lb;
              if (!silent)
              {
                  process_NE_repeat(process_NE_repeat_data,
                                  encseq,
                                  &suftab_arr.spaceGtUword[currelem.lb-
                                  (poscorr_flush+poscorr_mem)],
                                  currelem.lcp,
                                  ((idx+1)-(poscorr_flush+poscorr_mem))-
                                  (currelem.lb-(poscorr_flush+poscorr_mem)));
              }
            }  
          }
        }
      } 
      if (GT_STACK_TOP(&lcpstack).lcp < lcp && lcp >= searchlength)
      {
        currelem.lb = lb;
        currelem.lcp = lcp;
        GT_STACK_PUSH(&lcpstack,currelem);
      } 
    }

//    if (lcp < searchlength)
//    {
//      poscorr_mem--;
//    }

    if (lcp == 0)
    {
      poscorr_flush = idx;
      poscorr_mem = 0;
      suftab_arr.nextfreeGtUword = 0;
      GT_STOREINARRAY(&suftab_arr,GtUword,32,psuf);
    }
    plcp = lcp;
  }

  GT_FREEARRAY(&suftab_arr,GtUword);
  GT_STACK_DELETE(&lcpstack);
  if (nerepeat_progress != NULL)
  {
    gt_timer_show_progress_final(nerepeat_progress, stdout);
    gt_timer_delete(nerepeat_progress);
  }

  return haserr ? -1 : 0;
}

