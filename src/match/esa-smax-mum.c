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
#include "match/esa-smax-mum.h"

int gt_runlinsmaxmum(Sequentialsuffixarrayreader *ssar,
    GtUword searchlength,
    const bool silent,
    GtProcessSmaxpairs process_smaxpairs,
    void *process_smaxpairsdata,
    GtError *err)
{
  int haserr = 0;

  if (!haserr)
  {
    GtTimer *linsmaxprogress = NULL;
    const Suffixarray *suffixarray;
    GtUword lcpvalue,
              previoussuffix,
              idx,
              nonspecials,
              currentlcpmax,
              numseq;
    bool in_interval, prune;
    GtUchar bwt;
    GtBittab *marktab;
    GtBittab *mumtab;
    GtArrayGtUword suftab_arr;
    GtArrayGtUchar bwttab_arr;
    const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(ssar);

    gt_error_check(err);
    marktab = gt_bittab_new(gt_alphabet_size(gt_encseq_alphabet(encseq)));
    numseq = gt_encseq_num_of_sequences(encseq);
    mumtab = gt_bittab_new(numseq);
    nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
    suffixarray = gt_suffixarraySequentialsuffixarrayreader(ssar);
    GT_INITARRAY(&suftab_arr,GtUword);
    GT_INITARRAY(&bwttab_arr,GtUchar);

    gt_showtime_enable();
    if (gt_showtime_enabled())
    {
      linsmaxprogress = gt_timer_new_with_progress_description(
                                              "finding supermaximal repeats");
      gt_timer_start(linsmaxprogress);
    }
    bwt = 0;
    prune = true;
    in_interval = false;
    currentlcpmax = 0;
    for (idx = 0; idx < nonspecials; idx++)
    {
      SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
      SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
      gt_readnextfromstream_GtUchar(&bwt,(GtBufferedfile_GtUchar*)&suffixarray->bwttabstream);

      if (lcpvalue >= searchlength || !prune)
      {
        prune = false;

        if (lcpvalue > currentlcpmax)
        {
          in_interval = true;
          currentlcpmax = lcpvalue;
          suftab_arr.nextfreeGtUword = 0;
          bwttab_arr.nextfreeGtUchar = 0;
        }

        if (in_interval)
        {
          GT_STOREINARRAY(&suftab_arr,GtUword,32,previoussuffix);
          GT_STOREINARRAY(&bwttab_arr,GtUchar,32,bwt);
        }

        if (lcpvalue < currentlcpmax)
        {
          if (in_interval && suftab_arr.nextfreeGtUword <= numseq &&
                          gt_esa_smax_verify_mum(encseq,
                                      suftab_arr.spaceGtUword,
                                      suftab_arr.nextfreeGtUword,
                                      bwttab_arr.spaceGtUchar,
                                      bwttab_arr.nextfreeGtUchar,
                                      marktab,
                                      mumtab))
          {
            if (!silent)
            {
              process_smaxpairs (process_smaxpairsdata,
                                encseq,
                                currentlcpmax,
                                suftab_arr.spaceGtUword,
                                suftab_arr.nextfreeGtUword);
            }
          }          
          in_interval = false;
          currentlcpmax = lcpvalue;
          prune = true;
        }
      } else {
        currentlcpmax=lcpvalue;
      }
    }
    if (linsmaxprogress != NULL)
    {
      gt_timer_show_progress_final(linsmaxprogress, stdout);
      gt_timer_delete(linsmaxprogress);
    }
    GT_FREEARRAY(&suftab_arr,GtUword);
    GT_FREEARRAY(&bwttab_arr,GtUchar);
    gt_bittab_delete(marktab);
    gt_bittab_delete(mumtab);
  }
/*  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
    ssar = NULL;
  }
*/
  return haserr ? -1 : 0;
}