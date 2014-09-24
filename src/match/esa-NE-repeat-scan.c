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

/*
int gt_run_NE_repeats_scan(Sequentialsuffixarrayreader *ssar,
                          GT_UNUSED GtUword searchlength,
                          GT_UNUSED bool silent,
                          GT_UNUSED GtProcessNEintervals process_NE_repeat,
                          GT_UNUSED void *process_NE_repeat_data,
                          GtError *err)
{
  GT_UNUSED bool haserr = false;
  gt_error_check(err);
  GtTimer *nerepeat_progress = NULL;
  GtStackLcp_stackelem lcpstack;
//  Lcp_stackelem currelem;
  GtArrayGtUword suftab_arr;
  GtUword lcp,
          previoussuffix,
          idx,
          nonspecials,
          currentlcpmax;
  bool in_interval = true;
//  Definedunsignedint lctx_psuf,
//                     lctx_nsuf,
//                     lctx;
  GT_UNUSED const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(ssar);
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


  for (idx = 0; idx < nonspecials; idx++)
  {
    SSAR_NEXTSEQUENTIALLCPTABVALUE(lcp,ssar);
    SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    printf("SUF: " GT_WU "\t LCP: " GT_WU "\n",previoussuffix,lcp);

    if(lcpvalue > currentlcpmax)
    {
      in_interval = true;
      currentlcpmax = lcpvalue;
    }



  }
  return haserr ? -1 : 0;
}
*/

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
  Lcp_stackelem currelem;
  GtArrayGtUword suftab_arr;
  GtUword lcp,
          plcp,
          nsuf,
          psuf,
          idx,
          lb,
          poscorr_flush,
          poscorr_mem,
          nonspecials;
  Definedunsignedint lctx_psuf,
                     lctx_nsuf,
                     lctx;

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

//  lcp = 0;
  plcp = 0;
  psuf = 0;
//  lb = 0;
  poscorr_flush = 0;
  poscorr_mem = 0;

  for (idx = 0; idx < nonspecials; idx++)
  {
    if (idx == 0) 
    {
      SSAR_NEXTSEQUENTIALSUFTABVALUE(psuf,ssar);
      currelem.lcp = 0;
      currelem.lb = idx;
      if (psuf > 0)
      {
        lctx_psuf.defined = true;
        lctx_psuf = get_lctx(encseq,psuf);
        currelem.lctx = lctx_psuf;
      } else 
      {
        lctx_psuf.defined = false;
        lctx_psuf.valueunsignedint = 0;
        currelem.lctx = lctx_psuf;
      }
      GT_STOREINARRAY(&suftab_arr,GtUword,32,psuf);
      GT_STACK_PUSH(&lcpstack,currelem);
    } else
    {
      SSAR_NEXTSEQUENTIALLCPTABVALUE(lcp,ssar);
      SSAR_NEXTSEQUENTIALSUFTABVALUE(nsuf,ssar);
//printf("SUF: " GT_WU " LCP:" GT_WU "\n",nsuf,lcp);
      lb = idx-1;
      lctx_nsuf = get_lctx(encseq,nsuf);
      lctx = check_lctx(lctx_psuf, lctx_nsuf);
      lctx_psuf = lctx_nsuf;

      while (GT_STACK_TOP(&lcpstack).lcp > lcp)
      {
        currelem = GT_STACK_POP(&lcpstack);
        if (!currelem.lctx.defined && currelem.lcp >= searchlength)
        {
          if (!silent)
          {
            process_NE_repeat(process_NE_repeat_data,
                              encseq,
                              &suftab_arr.spaceGtUword[currelem.lb-
                              (poscorr_flush+poscorr_mem-1)],
                              currelem.lcp,
                              (idx-(poscorr_flush+poscorr_mem-1))-
                              (currelem.lb-(poscorr_flush+poscorr_mem-1)));
          }
        } 
        lb = currelem.lb;
        GT_STACK_TOP(&lcpstack).lctx = 
                      check_lctx(currelem.lctx,
                      GT_STACK_TOP(&lcpstack).lctx);
        lctx = check_lctx(currelem.lctx,lctx);
      }

      if (lcp >= searchlength) // decreases run-time
      {
        if (GT_STACK_TOP(&lcpstack).lcp == lcp)
        {
           GT_STOREINARRAY(&suftab_arr,GtUword,32,nsuf);
           GT_STACK_TOP(&lcpstack).lctx =
                                   check_lctx(
                                   GT_STACK_TOP(&lcpstack).lctx,
                                   lctx);
        } else 
        {
          if (plcp < searchlength)
          {
            poscorr_flush = idx;
            poscorr_mem = 0;
            suftab_arr.nextfreeGtUword = 0;
            GT_STOREINARRAY(&suftab_arr,GtUword,32,nsuf);
            GT_STOREINARRAY(&suftab_arr,GtUword,32,psuf);
            poscorr_mem--;
          }
          GT_STOREINARRAY(&suftab_arr,GtUword,32,nsuf); 
//          poscorr_mem++;
          currelem.lcp = lcp;
          currelem.lb = lb;
          currelem.lctx = lctx;
          GT_STACK_PUSH(&lcpstack,currelem);
        }
      } else {
        poscorr_mem++;
      }

      if (lcp == 0)
      {
        poscorr_flush = idx;
        poscorr_mem = 0;
        suftab_arr.nextfreeGtUword = 0;
        GT_STOREINARRAY(&suftab_arr,GtUword,32,nsuf);
      }
      plcp = lcp;
      psuf = nsuf;
    }
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

