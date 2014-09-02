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
//  GtUword suf;
  GtUword lb;
  Definedunsignedint lctx;
} Lcp_stackelem;

GT_STACK_DECLARESTRUCT(Lcp_stackelem, 32UL);

int gt_run_NE_repeats_map(Sequentialsuffixarrayreader *ssar,
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
  Lcp_stackelem currelem;
  const GtUword *suftab;
  GtUword lcp,
//          plcp,
          nsuf,
          psuf,
          idx,
          lb,
          nonspecials;
  Definedunsignedint lctx_psuf,
                     lctx_nsuf,
                     lctx;

  const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(ssar);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  suftab = gt_suftabSequentialsuffixarrayreader(ssar);
  GT_STACK_INIT(&lcpstack,32UL);

  gt_showtime_enable();
  if (gt_showtime_enabled())
  {
    nerepeat_progress = gt_timer_new_with_progress_description(
        "finding non-extendible repeats");
  gt_timer_start(nerepeat_progress);
  }

  lcp = 0;
//  plcp = 0;
//  nsuf = 0;
//  psuf = 0;
  lb = 0;

  for (idx = 0; idx < nonspecials; idx++)
  {
    if (idx == 0) 
    {
      SSAR_NEXTSEQUENTIALSUFTABVALUE(psuf,ssar);
      currelem.lcp = 0;
//      currelem.suf = psuf;
      currelem.lb = idx;
      if (psuf > 0)
      {
        lctx_psuf = get_lctx(encseq,psuf);
        currelem.lctx = lctx_psuf;
      } else 
      {
        lctx_psuf.defined = false;
        lctx_psuf.valueunsignedint = 0;
        currelem.lctx = lctx_psuf;
      }
      GT_STACK_PUSH(&lcpstack,currelem);
    } else
    {
      SSAR_NEXTSEQUENTIALLCPTABVALUE(lcp,ssar);
      SSAR_NEXTSEQUENTIALSUFTABVALUE(nsuf,ssar);
      lb = idx-1;

      lctx_nsuf = get_lctx(encseq,nsuf);
      lctx = check_lctx(lctx_psuf, lctx_nsuf);
      lctx_psuf = lctx_nsuf;

      while (GT_STACK_TOP(&lcpstack).lcp > lcp)
      {
        currelem = GT_STACK_POP(&lcpstack);
/*        printf("currelem.lctx: %c \n",
            currelem.lctx.valueunsignedint);
        printf("stack_top.lctx: %c \n",
                GT_STACK_TOP(&lcpstack).lctx.valueunsignedint);
*/
        if (!currelem.lctx.defined && currelem.lcp >= searchlength)
        {
          if (!silent)
          {
#ifdef SKDEBUG
            if(is_notleftextendible(encseq,
                                     &suftab[currelem.lb],
                                     idx-currelem.lb))
            {
              printf("Range is NLE\n");
            } else
            {
              printf("WARNING LE Range detected");
            }
            printf("LCP: " GT_WU " i: " GT_WU " j: " GT_WU "\n",
                  currelem.lcp, currelem.lb, idx-1);
#endif
            process_NE_repeat(process_NE_repeat_data,
                              encseq,
                              &suftab[currelem.lb],
                              currelem.lcp,
                              idx-currelem.lb);
          }
        } 
        lb = currelem.lb;
        GT_STACK_TOP(&lcpstack).lctx = 
                      check_lctx(currelem.lctx,
                        GT_STACK_TOP(&lcpstack).lctx);
        lctx = check_lctx(currelem.lctx,lctx);
      }

      if (GT_STACK_TOP(&lcpstack).lcp == lcp)
      {
            GT_STACK_TOP(&lcpstack).lctx =
                                    check_lctx(
                                    GT_STACK_TOP(&lcpstack).lctx,
                                    lctx);
      } else 
      {
        currelem.lcp = lcp;
//        currelem.suf = psuf;
        currelem.lb = lb;
        currelem.lctx = lctx;
        GT_STACK_PUSH(&lcpstack,currelem);
      }
/*
      printf("LCP: " GT_WU "\t",lcp);
      printf("PREVIOUSLCP: " GT_WU "\t",plcp);
      printf("LB: " GT_WU "\t",lb);
      printf("NSUF: " GT_WU "\t",nsuf);
      printf("PSUF: " GT_WU "\n",psuf);
*/
//      plcp = lcp;
//      psuf = nsuf;
    }
  }

  GT_STACK_DELETE(&lcpstack);
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
