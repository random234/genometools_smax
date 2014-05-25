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

/* vorraussichtlice Stackgroesse muss geschaetzt werden um
 * konstante Stackgroesse sinnvoll zu waehlen.
 */

typedef struct
{
  GtUword position;
  GtUword lcpvalue;
} GtLcp_stackelem;

GT_STACK_DECLARESTRUCT(GtLcp_stackelem, 32UL);

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
  GT_UNUSED GtLcp_stackelem lcpstack;
  GT_UNUSED GtUword lcpvalue,
            previoussuffix,
            idx,
            nonspecials,
            sumsuftab,
            sumlcptab,
            max = 0,
            period = -1;

  GT_UNUSED const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(ssar);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
//  GT_STACK_INIT(&lcpstack,32UL);

  gt_showtime_enable();
  if (gt_showtime_enabled())
  {
    nerepeat_progress = gt_timer_new_with_progress_description(
        "finding non-extendible repeats");
  gt_timer_start(nerepeat_progress);
  }

  for (idx = 0; idx < nonspecials; idx++)
  {
    SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
    sumlcptab += lcpvalue; /* silly but guarantees that loop is not
                              eliminated by compiler */
    SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    sumsuftab += previoussuffix; /* silly but guarantees that loop is not
                                  eliminated by compiler */

    while (period < max)
    {
      idx++;
      period = max;
      max = lcpvalue;
      if ((max > period) && (max >= searchlength))
      { 

        // push to stack
      }
    }

/*    while (top(LB).lcp <= max)
    {
      if (top(LB).lcp
    }
    */
  }

  if (nerepeat_progress != NULL)
  {
    gt_timer_show_progress_final(nerepeat_progress, stdout);
    gt_timer_delete(nerepeat_progress);
  }

  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
