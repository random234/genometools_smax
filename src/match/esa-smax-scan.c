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

static bool gt_linsmax_verify_supmax(Sequentialsuffixarrayreader *ssar,
    GtQueue *suftab_queue, GtESAMarkTab *marktab)
{
  GT_UNUSED const GtUword *suftab;
  GT_UNUSED const Suffixarray *sa;
  GtUword i;
  GtUword suf;
  sa = gt_suffixarraySequentialsuffixarrayreader(ssar);
  gt_esa_marktab_reset(marktab);

  /* printf("size of array: " GT_WU "\n", gt_queue_size(suftab_queue)); */
  /* push lb seq_pos
     * push lb seq_num */
/*  printf("verify_supmax called\n"); */
  for (i = 0;i<gt_queue_size(suftab_queue);i++)
  {
    suf = (GtUword) gt_queue_get(suftab_queue);
//    printf("SUF: " GT_WU "\t", suf); 
    if (suf > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(sa->encseq,suf-1,
          GT_READMODE_FORWARD);
//      printf("cc: %u \t",cc); 
      if (ISNOTSPECIAL(cc))
      {
        if (gt_esa_marktab_get(marktab,cc))
        {
          return false;
        }
        gt_esa_marktab_set(marktab,cc);
      }
    }
    gt_queue_add(suftab_queue, (GtUword*) suf);
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
//            nonspecials,
            seq_len,
            sumsuftab,
            sumlcptab,
            max;
  GtESAMarkTab *marktab;
  GtQueue *suftab_queue = gt_queue_new();

  sa = gt_suffixarraySequentialsuffixarrayreader(ssar);
  marktab = gt_esa_marktab_new(sa->encseq);
  seq_len = gt_Sequentialsuffixarrayreader_totallength(ssar);

/*  for (idx = 0; idx < seq_len; idx++)
  {
    SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
    sumlcptab += lcpvalue; 
    SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    sumsuftab += previoussuffix; 
    printf("SUF: \t" GT_WU "\t",previoussuffix);
    printf("LCP: \t" GT_WU "\n",lcpvalue);
  }
  
  ssar = NULL;
  ssar = gt_newSequentialsuffixarrayreaderfromfile(gt_str_array_get(
        inputindex,0),
        SARR_LCPTAB |
        SARR_SUFTAB |
        SARR_ESQTAB,
        false,
        logger,
        err);
*/
  GtUword queue_size;
  unsigned int i = 0;
  max = 0;
  for (idx = 0; idx < seq_len; idx++)
  {
    SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
    sumlcptab += lcpvalue; /* silly but guarantees that loop is not
                              eliminated by compiler */
    SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    sumsuftab += previoussuffix; /* silly but guarantees that loop is not
                                  eliminated by compiler */

      printf("current values: SUF: " GT_WU " LCP: " GT_WU "\n",
          previoussuffix, lcpvalue);
    if (lcpvalue > max)
    {
      /* flush queue */
      max = lcpvalue;
//      printf("max: " GT_WU "\n",max);
      queue_size = gt_queue_size(suftab_queue);      
      for (i = 0;i<queue_size;i++)
      { 
        gt_queue_get(suftab_queue);
      }
//      printf("queue flushed: " GT_WU "\n",gt_queue_size(suftab_queue));
    } 
    if (lcpvalue == max)
    {
      /* push to queue */
      gt_queue_add(suftab_queue, (GtUword*) previoussuffix);
//      printf("added value SUF: " GT_WU " to queue\n", previoussuffix);
    }
    if (lcpvalue < max)
    {
      /* push previoussuffix check local maxima and flush queue*/
//      printf("added value SUF: " GT_WU " to queue\n", previoussuffix);
      gt_queue_add(suftab_queue, (GtUword*) previoussuffix);
      if (gt_linsmax_verify_supmax(ssar, suftab_queue, marktab))
      {
//        printf("found supermaximal repeat\n");
        GtUword lb = (GtUword) gt_queue_get(suftab_queue);
        GtUword rb;
        queue_size = gt_queue_size(suftab_queue);
        for (i = 0;i<queue_size;i++)
        {
          rb = (GtUword) gt_queue_get(suftab_queue);
          printf("" GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " " GT_WU "\n",
                            max, lb, method, max, rb, 2*max);
          lb = rb;
        }
      } else
      {
//        printf("found local maxima, not supermaximal\n");
      }
//      printf("QUEUE SIZE: " GT_WU "\n",gt_queue_size(suftab_queue));
      queue_size = gt_queue_size(suftab_queue);
      for (i = 0;i<queue_size;i++)
      {
        gt_queue_get(suftab_queue);
//        printf("queue i: %d \n",i);
      }
//      printf("QUEUE flushed: " GT_WU "\n",gt_queue_size(suftab_queue));
      max = 0;
//      printf("max: " GT_WU "",max);
    }
  }

  if (linsmaxprogress != NULL)
  {
    /* gt_timer_show_progress(smaxprogress,"%GT_WD.%06GT_WDs real %GT_WDs user
     * %GT_WD system", stdout); */
    gt_timer_show_progress_final(linsmaxprogress, stdout);
    gt_timer_delete(linsmaxprogress);
  }

  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  gt_queue_delete(suftab_queue);
  gt_esa_marktab_delete(marktab);
  return haserr ? -1 : 0;
}
