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
#include "match/esa-seqread.h"
#include "match/esa-smax.h"
#include "match/esa-smax-scan.h"

int gt_runlinsmax(Sequentialsuffixarrayreader *ssar,
    GtUword searchlength,
    bool absolute,
    bool silent,
    GtProcessSmaxpairs process_smaxpairs,
    void *process_smaxpairsdata,
    GtError *err)
{
  bool haserr = false;
  
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
    
/*    printf("" GT_WU "\n",gt_encseq_num_of_sequences(encseq));
    for (idx = 0; idx <= nonspecials; idx++)
    {
        printf("%c",gt_encseq_get_decoded_char(encseq,idx,GT_READMODE_FORWARD));
    }
    printf("\n");
*/
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
          if (gt_esa_smax_verify_supmax(encseq, &suftab_arr, marktab))
          {
            GtUword s;
            for (s = 0;s<suftab_arr.nextfreeGtUlong;s++)
            {
              GtUword t;
              for (t = s+1;t<suftab_arr.nextfreeGtUlong;t++)
              {
                if (!silent)
                {
                  process_smaxpairs (process_smaxpairsdata,
                                      encseq,
                                      currentlcpmax,
                                      suftab_arr.spaceGtUlong[s],
                                      suftab_arr.spaceGtUlong[t],
                                      method,
                                      absolute);
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
