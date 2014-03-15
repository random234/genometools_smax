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

#include "lcpinterval.h"
#include "esa-smax-map.h"
#include "esa-seqread.h"
#include "esa-bottomup.h"
#include "core/encseq_api.h"
#include "core/str_api.h"
#include "esa-dfs.h"
#include "core/showtime.h"

int gt_runsmaxlcpvalues(GtStrArray *inputindex,
    GtUword searchlength,
    bool absolute,
    bool silent,
    GT_UNUSED bool outedges,
    bool bottomup,
    GtLogger *logger,
    GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;
  GtTimer *smaxprogress = NULL;
  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(gt_str_array_get(
        inputindex,0),
      SARR_LCPTAB |
      SARR_SUFTAB |
      SARR_ESQTAB,
      false,
      /*SEQ_scan,*/
      logger,
      err);

  gt_showtime_enable();
  if (gt_showtime_enabled())
  {
    smaxprogress = gt_timer_new_with_progress_description(
        "finding supermaximal repeats");
    gt_timer_start(smaxprogress);
  }

  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (bottomup)
    {
      GtESAVisitor *elv = gt_esa_smax_lcpitvs_visitor_new(
          ssar, searchlength, absolute, silent, smaxprogress);
      if (gt_esa_bottomup(ssar, elv, err) != 0)
      {
        haserr = true;
      }
      gt_esa_visitor_delete(elv);
    }
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }

  if (smaxprogress != NULL)
  {
/* gt_timer_show_progress(smaxprogress,"%GT_WD.%06GT_WDs real %GT_WDs
 * user %GT_WDs system", stdout); */
    gt_timer_show_progress_final(smaxprogress, stdout);
    gt_timer_delete(smaxprogress);
  }
  return haserr ? -1 : 0;
}
