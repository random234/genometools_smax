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
gt_log_log("Leaf exit_code: %c fd: GT_WU flb: GT_WU leafnumber: GT_WU\n",
    firstsucc ? '1' : '0', fd, flb, leafnumber);
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
  GtUword pos_corr_t, pos_corr_s;
  GtLcpmaxintervalinfo *ret_info = (GtLcpmaxintervalinfo *) info;

  gt_assert(ev && err);
  lev = gt_esa_smax_lcpitvs_visitor_cast(ev);
  sa = gt_suffixarraySequentialsuffixarrayreader(lev->ssar);
  suftab = gt_suftabSequentialsuffixarrayreader(lev->ssar);

  if (lcp >= lev->searchlength && ret_info->maxlcpinterval)
  {
    if (verify_supmax(lev,lb,rb)) {
      for (s=lb;s<rb;s++)
      {
        GtUword seqnum_s = gt_encseq_seqnum(sa->encseq,suftab[s]);
        for (t=s+1;t<=rb;t++)
        {
          if (suftab[s] < suftab[t]) 
            /* Da wir nicht wissen wer groesser oder kleiner ist existieren hier di Fallunterscheidung die dafuer sorgt das in jedem Fall die richtige Ausrichtung der suftab Tabelle genutzt wird */
          {
            /* Laenge des Repeats, Startpos, Methode(Direct), Laenge des korrespondierenden Repeats, Startpos 
              Sequenznummer wird nur bei relativer Positionsangabe verwendet*/
            if (lev->absolute)
            {
              if (!lev->silent)
              {
                printf("%2lu %3lu %3c %2lu %3lu %3lu\n",lcp, suftab[s], method, lcp, suftab[t],lcp+lcp);
              }
            } else
            {
// return start position of seqnumth position: gt_encseq_seqstartpos(const GtEncseq *encseq, GtUword seqnum);
// printf("Sequence number: GT_WU \n", gt_encseq_seqnum(sa->encseq,suftab[t]));
// printf("Sequence Start position: GT_WU \n", gt_encseq_seqstartpos(sa->encseq, gt_encseq_seqnum(sa->encseq,suftab[t])));
              pos_corr_t = gt_encseq_seqstartpos(sa->encseq, gt_encseq_seqnum(sa->encseq,suftab[t]));
              pos_corr_s = gt_encseq_seqstartpos(sa->encseq, seqnum_s);
              if (!lev->silent)
              {
                printf("%2lu %4lu %3lu %3c %2lu %4lu %3lu %3lu\n",lcp, gt_encseq_seqnum(sa->encseq,suftab[s]), suftab[s]-pos_corr_s, method, lcp, gt_encseq_seqnum(sa->encseq,suftab[t]), suftab[t]-pos_corr_t,lcp+lcp);
              }
            }
          } else
          {
            if (lev->absolute)
            {
              if (!lev->silent)
              {
                printf("%2lu %3lu %3c %2lu %3lu %3lu\n",lcp, suftab[t], method, lcp, suftab[s],lcp+lcp);
              }
            }
            else
            {
// printf("Sequence number: GT_WU \n", gt_encseq_seqnum(sa->encseq,suftab[t]));
// printf("Sequence Start position: GT_WU \n", gt_encseq_seqstartpos(sa->encseq, gt_encseq_seqnum(sa->encseq,suftab[t])));
              pos_corr_t = gt_encseq_seqstartpos(sa->encseq, gt_encseq_seqnum(sa->encseq,suftab[t]));
              pos_corr_s = gt_encseq_seqstartpos(sa->encseq, gt_encseq_seqnum(sa->encseq,suftab[s]));
              if(!lev->silent)
              {
                printf("%2lu %4lu %3lu %3c %2lu %4lu %3lu %3lu\n",lcp, gt_encseq_seqnum(sa->encseq,suftab[t]), suftab[t]-pos_corr_t, method, lcp, gt_encseq_seqnum(sa->encseq,suftab[s]), suftab[s]-pos_corr_s,lcp+lcp);
              }
            }
          }
        }
      }
    }
  }
  return 0;
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
    if (suftab[i] > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(sa->encseq,suftab[i]-1,
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
