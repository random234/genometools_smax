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
#include "core/bittab_api.h"
#include "core/arraydef.h"
#include "esa-seqread.h"
#include "esa-smax-callbacks.h"
#include "match/esa-smax-scan.h"

static bool gt_smaxscan_verify_supmax(const GtEncseq *encseq,
    const GtArrayGtUlong *suftab_arr, GtBittab *marktab)
{
  GtUword i;
  gt_bittab_unset(marktab);
  
  for (i = 0;i<suftab_arr->nextfreeGtUlong;i++)
  {
    GtUword suf = suftab_arr->spaceGtUlong[i];
    if (suf > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,suf-1,
                                              GT_READMODE_FORWARD);
      if (ISNOTSPECIAL(cc))
      {
        if (gt_bittab_bit_is_set(marktab,cc))
        {
          return false;
        }
        gt_bittab_set_bit(marktab,cc);
      }
    }
  }
  return true;
}

/*
void gt_smaxscan_print_repeat(const GtEncseq *encseq, GtUword maxlen,
                              GtUword suftab_s, GtUword suftab_t,
                              char method, bool absolute)
{
  GtUword score = maxlen * 2;
  GtUword seqnum_s = gt_encseq_seqnum(encseq,suftab_s);
  GtUword seqnum_t = gt_encseq_seqnum(encseq,suftab_t);
  if (suftab_s > suftab_t)
  {
    GtUword tmp;
    tmp = suftab_s;
    suftab_s = suftab_t;
    suftab_t = tmp;
    tmp = seqnum_s;
    seqnum_s = seqnum_t;
    seqnum_t = tmp;
  }
  if (absolute)
  {
    printf(""GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " " GT_WU"\n",maxlen,
        suftab_s, method, maxlen, suftab_t,score);
  } else
  {
    GtUword pos_corr_t = gt_encseq_seqstartpos(encseq, seqnum_t),
            pos_corr_s = gt_encseq_seqstartpos(encseq, seqnum_s);

    printf("" GT_WU " " GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " "
        GT_WU " " GT_WU "\n",maxlen, seqnum_s, suftab_s-pos_corr_s, method,
        maxlen, seqnum_t, suftab_t-pos_corr_t,score);
  }
}
*/


int gt_runlinsmax(GtStrArray *inputindex,
    GtUword searchlength,
    bool absolute,
    bool silent,
    GT_UNUSED bool outedges,
    GtLogger *logger,
    GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar =
    gt_newSequentialsuffixarrayreaderfromfile(gt_str_array_get(
        inputindex,0),
        SARR_LCPTAB |
        SARR_SUFTAB |
        SARR_ESQTAB,
        true, /* scan suftab and lcptab */
        logger,
        err);

  if (ssar == NULL)
  {
    haserr = true;
  }  
  if (!haserr)
  {
    GtTimer *linsmaxprogress = NULL;
    char method = 'D';
    GtUword lcpvalue,
              previoussuffix,
              idx,
              nonspecials,
              currentlcpmax;
    GtBittab *marktab;
    GtArrayGtUlong suftab_arr;
    const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(ssar);
    bool in_interval = true;
    gt_error_check(err);
    marktab = gt_bittab_new(gt_alphabet_size(gt_encseq_alphabet(encseq)));
    nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
    GT_INITARRAY(&suftab_arr,GtUlong);

    gt_showtime_enable();
    if (gt_showtime_enabled())
    {
      linsmaxprogress = gt_timer_new_with_progress_description(
                                              "finding supermaximal repeats");
      gt_timer_start(linsmaxprogress);
    }

    currentlcpmax = 0;
    for (idx = 0; idx < nonspecials; idx++)
    {
      SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
      SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
      if (lcpvalue > currentlcpmax)
      {
        in_interval = true;
        currentlcpmax = lcpvalue;
        suftab_arr.nextfreeGtUlong = 0;
      }
      if (lcpvalue == currentlcpmax && in_interval)
      {
        GT_STOREINARRAY(&suftab_arr,GtUlong,32,previoussuffix);
      }
      if (lcpvalue < currentlcpmax)
      {
        if (in_interval && currentlcpmax >= searchlength)
        {
          GT_STOREINARRAY(&suftab_arr,GtUlong,32,previoussuffix);
          if (gt_smaxscan_verify_supmax(encseq, &suftab_arr, marktab))
          {
            GtUword s;
            for (s = 0;s<suftab_arr.nextfreeGtUlong;s++)
            {
              GtUword t;
              for (t = s+1;t<suftab_arr.nextfreeGtUlong;t++)
              {
                if (!silent)
                {
                  gt_esa_smax_print_repeat(encseq, currentlcpmax,
                                          suftab_arr.spaceGtUlong[s],
                                          suftab_arr.spaceGtUlong[t],
                                          method, absolute);
                }
              }
            }
          }
        }
        in_interval = false;
        currentlcpmax = lcpvalue;
        suftab_arr.nextfreeGtUlong = 0;
      }
    }

    if (linsmaxprogress != NULL)
    {
      gt_timer_show_progress_final(linsmaxprogress, stdout);
      gt_timer_delete(linsmaxprogress);
    }
    GT_FREEARRAY(&suftab_arr,GtUlong);
    gt_bittab_delete(marktab);
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
