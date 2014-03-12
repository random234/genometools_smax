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

#include "core/class_alloc_lock.h"
#include "core/unused_api.h"
#include "core/log_api.h"
#include "esa_visitor_rep.h"
#include "esa-smax-map.h"
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
    GtESAVisitorInfo *info,
    GT_UNUSED GtUword leafnumber,
    GT_UNUSED GtError *err)
{
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
    GtESAVisitorInfo *finfo,
    GT_UNUSED GtUword sd,
    GT_UNUSED GtUword slb,
    GT_UNUSED GtUword srb,
    GT_UNUSED GtESAVisitorInfo *sinfo,
    GT_UNUSED GtError *err)
{
  GtLcpmaxintervalinfo *ret_info = (GtLcpmaxintervalinfo *) finfo;
  ret_info->maxlcpinterval = false;
  return 0;
}

static bool verify_supmax(GtESASmaxLcpintervalsVisitor *lev, GtUword lb,
                          GtUword rb)
{
  GtUword idx;
  const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(lev->ssar);
  const GtUword *suftab = gt_suftabSequentialsuffixarrayreader(lev->ssar);

  gt_esa_marktab_reset(lev->marktab);
  for (idx=lb;idx<=rb;idx++)
  {
    if (suftab[idx] > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,suftab[idx]-1,
                                              GT_READMODE_FORWARD);
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

void gt_esa_smax_print_repeat(const GtEncseq *encseq, GtUword maxlen,
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

static int gt_esa_smax_lcpitvs_visitor_visitlcpinterval(GtESAVisitor *ev,
    GtUword lcp,
    GtUword lb,
    GtUword rb,
    GtESAVisitorInfo *info,
    GtError *err)
{
  GtESASmaxLcpintervalsVisitor *lev;
  GtUword s,t;
  char method = 'D';
  GtLcpmaxintervalinfo *ret_info = (GtLcpmaxintervalinfo *) info;
  gt_assert(ev && err);

  lev = gt_esa_smax_lcpitvs_visitor_cast(ev);
  const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(lev->ssar);
  const GtUword *suftab = gt_suftabSequentialsuffixarrayreader(lev->ssar);

  if (lcp >= lev->searchlength && ret_info->maxlcpinterval)
  {
    if (verify_supmax(lev,lb,rb))
    {
      for (s=lb;s<rb;s++)
      {
        for (t=s+1;t<=rb;t++)
        {
          if (!lev->silent)
          {
            gt_esa_smax_print_repeat(encseq, lcp, suftab[s], suftab[t],
                method, lev->absolute);
          }
        }
      }
    }
  }
  return 0;
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
