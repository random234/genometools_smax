/*
  Copyright (c) 2011 Sascha Steinbiss <random234@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/class_alloc_lock.h"
#include "core/unused_api.h"
#include "core/log_api.h"
#include "esa_visitor_rep.h"
#include "esa_smax_lcpvisitor.h"
#include "esa-seqread.h"
#include "esa-marktab.h"

struct GtESASmaxLcpintervalsVisitor {
  const GtESAVisitor parent_instance;
  Sequentialsuffixarrayreader *ssar;
  GtUword searchlength;
  GtTimer *smaxprogress;
  bool absolute;
  bool silent;
  GtESAMarkTab *marktab;
};

typedef struct
{
  bool maxlcpinterval;
} GtLcpmaxintervalinfo;

#define gt_esa_smax_lcpitvs_visitor_cast(GV)\
        gt_esa_visitor_cast(gt_esa_smax_lcpitvs_visitor_class(), GV)

static GtESAVisitorInfo *gt_esa_smax_lcpitvs_visitor_create_info(
    GT_UNUSED GtESAVisitor *ev)
{
  GtLcpmaxintervalinfo *info = gt_malloc(sizeof(*info));
  info->maxlcpinterval = true;
  return (GtESAVisitorInfo *) info;
}

static void gt_esa_smax_lcpitvs_visitor_delete_info(GtESAVisitorInfo *vi,
                                                GtESAVisitor *ev)
{
  GtESASmaxLcpintervalsVisitor *lev = gt_esa_smax_lcpitvs_visitor_cast(ev);
  gt_esa_marktab_delete(lev->marktab);
  lev->marktab = NULL;
  gt_free(vi);
}

static int gt_esa_smax_lcpitvs_visitor_processleafedge(
    GT_UNUSED GtESAVisitor *ev,
    GT_UNUSED bool firstsucc,
    GT_UNUSED GtUword fd,
    GT_UNUSED GtUword flb,
    GT_UNUSED GtESAVisitorInfo *info,
    GT_UNUSED GtUword leafnumber,
    GT_UNUSED GtError *err)
{
gt_log_log("Leaf exit_code: %c fd:" GT_WU " flb:" GT_WU " leafnumber:" GT_WU
    "\n", firstsucc ? '1' : '0', fd, flb, leafnumber);
  if (firstsucc)
  {
    ((GtLcpmaxintervalinfo *) info)->maxlcpinterval = true;
  }
  return 0;
}
static int gt_esa_smax_lcpitvs_visitor_processbranchingedge(
    GT_UNUSED GtESAVisitor *ev,
    GT_UNUSED bool firstsucc,
    GT_UNUSED GtUword fd,
    GT_UNUSED GtUword flb,
    GT_UNUSED GtESAVisitorInfo *finfo,
    GT_UNUSED GtUword sd,
    GT_UNUSED GtUword slb,
    GT_UNUSED GtUword srb,
    GT_UNUSED GtESAVisitorInfo *sinfo,
    GT_UNUSED GtError *err)
{
  GtLcpmaxintervalinfo *ret_info = (GtLcpmaxintervalinfo *) finfo;
  ret_info->maxlcpinterval = false;
//  gt_timer_show_progress(lev->smaxprogress,"visiting branch", stdout);
//  gt_esa_smax_lcpitvs_visitor_set_info(finfo, false);
  return 0;
}

static int gt_esa_smax_lcpitvs_visitor_visitlcpinterval(GtESAVisitor *ev,
    GtUword lcp,
    GtUword lb,
    GtUword rb,
    GtESAVisitorInfo *info,
    GtError *err)
{
  GtESASmaxLcpintervalsVisitor *lev;
  const GtUword *suftab;
  const Suffixarray *sa;
  GtUword s,t;
  char method = 'D';
  GtLcpmaxintervalinfo *ret_info = (GtLcpmaxintervalinfo *) info;

  gt_assert(ev && err);
  lev = gt_esa_smax_lcpitvs_visitor_cast(ev);
  sa = gt_suffixarraySequentialsuffixarrayreader(lev->ssar);
  suftab = gt_suftabSequentialsuffixarrayreader(lev->ssar);

  if (lcp >= lev->searchlength && ret_info->maxlcpinterval)
  {
    if (verify_supmax(lev,lb,rb))
    {
      for (s=lb;s<rb;s++)
      {
        GtUword seqnum_s = gt_encseq_seqnum(sa->encseq,suftab[s]);
        for (t=s+1;t<=rb;t++)
        {
          GtUword seqnum_t = gt_encseq_seqnum(sa->encseq,suftab[t]);
          if (!lev->silent)
          {
            gt_esa_smax_print_repeat(sa->encseq, lcp, suftab[s], suftab[t],
                seqnum_s, seqnum_t, method, lev->absolute);
          }
        }
      }
    }
  }
  return 0;
}

void gt_esa_smax_print_repeat(GtEncseq *encseq, GtUword lcp, GtUword suftab_s,
    GtUword suftab_t, GtUword seqnum_s, GtUword seqnum_t,  char method,
    bool absolute)
{
  GtUword pos_corr_s, pos_corr_t;
  if (suftab_s > suftab_t)
  {
    gt_esa_smax_print_repeat(encseq, lcp,suftab_t,suftab_s, seqnum_t, seqnum_s,
        method, absolute);
  }

  if (absolute)
  {
    printf(""GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " " GT_WU"\n",lcp,
        suftab_s, method, lcp, suftab_t,lcp+lcp);
  } else
  {
    pos_corr_t = gt_encseq_seqstartpos(encseq, seqnum_t);
    pos_corr_s = gt_encseq_seqstartpos(encseq, seqnum_s);

    printf("" GT_WU " " GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " " GT_WU " "
        GT_WU "\n",lcp, seqnum_s, suftab_s-pos_corr_s, method, lcp, seqnum_t,
        suftab_t-pos_corr_t,lcp+lcp);
  }
}

bool verify_supmax(GtESASmaxLcpintervalsVisitor *lev, GtUword lb, GtUword rb) {
  const GtUword *suftab;
  const Suffixarray *sa;
/*  bool marktab[GT_DNAALPHASIZE]; nur einmal allokieren und zwar
 *  entsprechend der Alphabetgr"osse */
  int i;
  sa = gt_suffixarraySequentialsuffixarrayreader(lev->ssar);
  suftab = gt_suftabSequentialsuffixarrayreader(lev->ssar);
  gt_esa_marktab_reset(lev->marktab);
  for (i=lb;i<=rb;i++)
  {
//    printf("SUF; " GT_WU "\n",suftab[i]);
    if (suftab[i] > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(sa->encseq,suftab[i]-1,
          GT_READMODE_FORWARD);
//      printf("cc: %u\n",cc);
      if (ISNOTSPECIAL(cc))
      {
        if (gt_esa_marktab_get(lev->marktab,cc))
        {
          return false;
        }
        gt_esa_marktab_set(lev->marktab,cc);
      }
    }
  }
  return true;
}

const GtESAVisitorClass* gt_esa_smax_lcpitvs_visitor_class()
{
  static const GtESAVisitorClass *esc = NULL;
  gt_class_alloc_lock_enter();
  if (!esc) {
    esc = gt_esa_visitor_class_new(sizeof (GtESASmaxLcpintervalsVisitor),
    NULL,
    gt_esa_smax_lcpitvs_visitor_processleafedge,
    gt_esa_smax_lcpitvs_visitor_processbranchingedge,
    gt_esa_smax_lcpitvs_visitor_visitlcpinterval,
    gt_esa_smax_lcpitvs_visitor_create_info,
    gt_esa_smax_lcpitvs_visitor_delete_info);
  }
  gt_class_alloc_lock_leave();
  return esc;
}

GtESAVisitor* gt_esa_smax_lcpitvs_visitor_new(
    Sequentialsuffixarrayreader *ssar, GtUword searchlength, bool absolute,
    bool silent, GtTimer *smaxprogress)
{
  GtESASmaxLcpintervalsVisitor *lev;
  GtESAVisitor *ev = gt_esa_visitor_create(gt_esa_smax_lcpitvs_visitor_class());
  lev = gt_esa_smax_lcpitvs_visitor_cast(ev);
  lev->ssar = ssar;
  lev->searchlength = searchlength;
  lev->absolute = absolute;
  lev->smaxprogress = smaxprogress;
  lev->silent = silent;
  lev->marktab = gt_esa_marktab_new(ssar->encseq);
  return ev;
}
