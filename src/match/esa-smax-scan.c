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
#include "esa-smax.h"
#include "match/esa-smax-scan.h"

typedef struct {
  const GtEncseq *encseq;
  GtArrayGtUlong *suftab_arr;
  GtBittab *marktab;
} GtESASmaxScanVerifyInput;


static bool scan_verify_supmax(void *data)
{
  GtESASmaxScanVerifyInput *input = (GtESASmaxScanVerifyInput*) data; 
  GtUword i;
  gt_bittab_unset(input->marktab);

  for (i = 0;i<input->suftab_arr->nextfreeGtUlong;i++)
  {
    GtUword suf = input->suftab_arr->spaceGtUlong[i];
    if (suf > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(input->encseq,suf-1,
                                             GT_READMODE_FORWARD);
      if (ISNOTSPECIAL(cc))
      {
        if (gt_bittab_bit_is_set(input->marktab,cc))
        {
          return false;
        }
        gt_bittab_set_bit(input->marktab,cc);
      }
    }
  }
  return true;
}

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
    GtESASmaxScanVerifyInput *input = gt_malloc(sizeof(*input));
    input->encseq = encseq;

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
          input->suftab_arr = &suftab_arr;
          input->marktab = marktab;
          if (gt_esa_smax_verify_supmax(scan_verify_supmax,input))
          /*if (gt_smaxscan_verify_supmax(encseq, &suftab_arr, marktab))*/
          {
            GtUword s;
            for (s = 0;s<suftab_arr.nextfreeGtUlong;s++)
            {
              GtUword t;
              for (t = s+1;t<suftab_arr.nextfreeGtUlong;t++)
              {
                if (!silent)
                {
                  print_repeat_both (NULL, encseq,
                                    currentlcpmax, 
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
    gt_free(input);
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
