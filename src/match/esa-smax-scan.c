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

#include "esa-seqread.h"
#include "core/encseq_api.h"
#include "core/str_api.h"
#include "core/showtime.h"
#include "match/esa-smax-scan.h"
#include "core/queue_api.h"
#include "esa-marktab.h"

static bool gt_linsmax_verify_supmax(const GtEncseq *encseq,
    GtArray *suftab_arr, GtESAMarkTab *marktab)
{
  GtUword i;
  GtUword suf;
  gt_esa_marktab_reset(marktab);

  for (i = 0;i<gt_array_size(suftab_arr);i++)
  {
    suf = *(GtUword*) gt_array_get(suftab_arr,i);
    if (suf > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,suf-1,
          GT_READMODE_FORWARD);
      if (ISNOTSPECIAL(cc))
      {
        if (gt_esa_marktab_get(marktab,cc))
        {
          return false;
        }
        gt_esa_marktab_set(marktab,cc);
      }
    }
  }
  return true;
}

int gt_runlinsmax(GtStrArray *inputindex,
    GT_UNUSED GtUword searchlength,
    GT_UNUSED bool absolute,
    GT_UNUSED bool silent,
    GT_UNUSED bool outedges,
    GtLogger *logger,
    GtError *err)
{
  bool haserr = false;
  gt_error_check(err);
  Sequentialsuffixarrayreader *ssar =
    gt_newSequentialsuffixarrayreaderfromfile(gt_str_array_get(
        inputindex,0),
        SARR_LCPTAB |
        SARR_SUFTAB |
        SARR_ESQTAB,
        false,
        /* SEQ_scan, */
        logger,
        err);
  GtTimer *linsmaxprogress = NULL;
  GT_UNUSED const Suffixarray *sa;
  GT_UNUSED char method = 'D';
  GT_UNUSED int pos_corr_t,pos_corr_s = 0;
  GT_UNUSED GtUword lcpvalue,
            previoussuffix,
            idx,
            nonspecials,
            seq_len,
            sumsuftab,
            sumlcptab,
            max;
  GtESAMarkTab *marktab;
  GtArray *suftab_arr = gt_array_new(sizeof(GtUword));
  bool in_interval = true;

  const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(ssar);
  marktab = gt_esa_marktab_new(encseq);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);

  gt_showtime_enable();
  if (gt_showtime_enabled())
  {
    linsmaxprogress = gt_timer_new_with_progress_description(
        "finding supermaximal repeats");
  gt_timer_start(linsmaxprogress);
  }

  unsigned int s,t;
  max = 0;
  for (idx = 0; idx < nonspecials; idx++)
  {
    SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
    sumlcptab += lcpvalue; /* silly but guarantees that loop is not
                              eliminated by compiler */
    SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    sumsuftab += previoussuffix; /* silly but guarantees that loop is not
                                  eliminated by compiler */
    /* check if we are currently in a lcp-interval */
    if (lcpvalue > max)
    {
      /* flush queue */
      in_interval = true;
      max = lcpvalue;
      gt_array_reset(suftab_arr);
    } 
    if (lcpvalue == max && in_interval)
    {
      /* push to queue */
      gt_array_add(suftab_arr,previoussuffix);
    }
    if (lcpvalue < max)
    {
      if (in_interval)
      {
        /* push previoussuffix check local maxima and flush queue*/
        gt_array_add(suftab_arr,previoussuffix);
        if (gt_linsmax_verify_supmax(encseq, suftab_arr, marktab))
        {
          for (s = 0;s<gt_array_size(suftab_arr);s++)          
          {
            for (t = s+1;t<gt_array_size(suftab_arr);t++)
            {
             printf("" GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " " GT_WU "\n",
                  max, *(GtUword*) gt_array_get(suftab_arr,s), method, max,
                  *(GtUword*) gt_array_get(suftab_arr,t), max*2);
            }
          }
        }
      }
      if (max > lcpvalue)
      {
        in_interval = false;
      }
      gt_array_reset(suftab_arr);
      max = lcpvalue;
    }
  }

  if (linsmaxprogress != NULL)
  {
/*     gt_timer_show_progress(smaxprogress,"%GT_WD.%06GT_WDs real %GT_WDs user
     * %GT_WD system", stdout); */
    gt_timer_show_progress_final(linsmaxprogress, stdout);
    gt_timer_delete(linsmaxprogress);
  }

  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  gt_array_delete(suftab_arr);
  gt_esa_marktab_delete(marktab);
  return haserr ? -1 : 0;
}
