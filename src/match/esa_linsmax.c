#include "esa-seqread.h"
#include "core/encseq_api.h"
#include "core/str_api.h"
#include "core/showtime.h"
#include "match/esa_linsmax.h"
#include "core/queue_api.h"

int gt_runlinsmax(GtStrArray *inputindex,
    GT_UNUSED GtUword searchlength,
    GT_UNUSED bool absolute,
    GT_UNUSED bool silent,
    GT_UNUSED bool outedges,
    GtLogger *logger,
    GtError *err)
{
  gt_assert(inputindex);
  bool haserr = false;
  GT_UNUSED Sequentialsuffixarrayreader *ssar;
  GtTimer *linsmaxprogress = NULL;
  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(gt_str_array_get(
        inputindex,0),
        SARR_LCPTAB |
        SARR_SUFTAB |
        SARR_ESQTAB,
        false,
        /* SEQ_scan, */
        logger,
        err);
  GT_UNUSED const Suffixarray *sa;
  char method = 'D';
  GT_UNUSED int pos_corr_t,pos_corr_s = 0;
  GT_UNUSED GtUword lcpvalue,
            previoussuffix,
            idx,
            nonspecials,
            sumsuftab,
            sumlcptab,
            max,
            lb,
            lb_seq_pos,
            lb_seq_num,
            rb,
            rb_seq_pos,
            rb_seq_num,
            i;

  GtQueue *suftab_queue = NULL;
  suftab_queue = gt_queue_new();

  sa = gt_suffixarraySequentialsuffixarrayreader(ssar);

  gt_showtime_enable();
  if (gt_showtime_enabled())
  {
    linsmaxprogress = gt_timer_new_with_progress_description
      ("finding supermaximal repeats");
   gt_timer_start(linsmaxprogress);
  }
  
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  max = 0;
  for (idx = 0; idx < nonspecials; idx++)
  {
    SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
    sumlcptab += lcpvalue; /* silly but guarantees that loop is not eliminated
                            by compiler */
    SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    sumsuftab += previoussuffix; /* silly but guarantees that loop is not
                                  eliminated by compiler */

    gt_queue_add(suftab_queue, (GtUword*) previoussuffix);
    if (lcpvalue > max)
    {
      max = lcpvalue;
/* SOLANGE WEITERE MAXIMALE ELEMENTE GEFUNDEN WERDEN ELIMINIERE DIE
 * VORGAENGER */
      if (gt_queue_size(suftab_queue) > 2)
      {
        gt_queue_get(suftab_queue);
      }
    } else
    {
      /* gt_queue_add(suftab_queue, (GtUword*) previoussuffix); */
      if (max>=searchlength)
      {
        if (gt_linsmax_verify_supmax(ssar, suftab_queue, lb, rb))
        {
        /* Laenge des Repeats, Startpos, Methode(Direct), Laenge des
         * korrespondierenden Repeats, Startpos Sequenznummer wird nur bei
         * relativer Positionsangabe verwendet */
          lb = (GtUword) gt_queue_get(suftab_queue);
          rb = (GtUword) gt_queue_get(suftab_queue);
          printf("" GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " " GT_WU "\n",max, lb, method,
              max, rb, max+max);
          max=0;

          for (i = 0;i<gt_queue_size(suftab_queue);i++)
          {
            gt_queue_get(suftab_queue);
          }
        } else
        {
          /* printf("found local maxima, not supermaximal\n"); */
          max=0;
          for (i = 0;i<gt_queue_size(suftab_queue);i++)
          {
            gt_queue_get(suftab_queue);
          }
        }
      }
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
  return haserr ? -1 : 0;
}

bool gt_linsmax_verify_supmax(Sequentialsuffixarrayreader *ssar,
    GtQueue *suftab_queue, GT_UNUSED GtUword lb, GT_UNUSED GtUword rb)
{
  GT_UNUSED const GtUword *suftab;
  GT_UNUSED const Suffixarray *sa;
  GT_UNUSED bool marktab[GT_DNAALPHASIZE+1];
  GtUword i;
  GtUword suf;

  sa = gt_suffixarraySequentialsuffixarrayreader(ssar);
  /* suftab = gt_suftabSequentialsuffixarrayreader(ssar); */

  for (i=0;i<=GT_DNAALPHASIZE;i++)
  {
    marktab[i] = false;
  }

  /* printf("size of array: " GT_WU "\n", gt_queue_size(suftab_queue)); */
  for (i = 0;i<gt_queue_size(suftab_queue);i++) {
    suf = (GtUword) gt_queue_get(suftab_queue);
    gt_queue_add(suftab_queue, (GtUword*) suf);
    /* push lb seq_pos
     * push lb seq_num */
  }

  for (i = 0;i<gt_queue_size(suftab_queue);i++)
  {
    suf = (GtUword) gt_queue_get(suftab_queue);
    if (suf!=0)
    {
      if (ISNOTSPECIAL(gt_encseq_get_encoded_char(sa->encseq,suf-1,0)))
      {
        if (marktab[gt_encseq_get_encoded_char(sa->encseq,suf-1,0)] == true)
        {
          return false;
        }
        marktab[gt_encseq_get_encoded_char(sa->encseq,suf-1,0)] = true;
      }
    }
    gt_queue_add(suftab_queue, (GtUword*) suf);
  }
  return true;
}
