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
#include "esa_visitor_rep.h"
#include "esa-seqread.h"
#include "core/bittab_api.h"
#include "esa-smax.h"
#include "esa-smax-map.h"

typedef struct {
  GtESASmaxLcpintervalsVisitor *visitor;
  GtUword lb;
  GtUword rb;
} GtESASmaxMapVerifyInput;

struct GtESASmaxLcpintervalsVisitor {
  const GtESAVisitor parent_instance;
  Sequentialsuffixarrayreader *ssar;
  GtUword searchlength;
  GtTimer *smaxprogress;
  bool absolute;
  bool silent;
  GtBittab *marktab;
  GtESASmaxMapVerifyInput *input;
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
  gt_bittab_delete(lev->marktab);
  lev->marktab = NULL;
  gt_free(lev->input);
  lev->input = NULL;
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

static bool map_verify_supmax(void *data)
/*    GtESASmaxLcpintervalsVisitor *lev, GtUword lb,GtUword rb) */
{
  GtESASmaxMapVerifyInput *input = (GtESASmaxMapVerifyInput*) data;

  GtUword idx;
  const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(
                                                      input->visitor->ssar);
  const GtUword *suftab = gt_suftabSequentialsuffixarrayreader(
                                                      input->visitor->ssar);

  gt_bittab_unset(input->visitor->marktab);
  for (idx=input->lb;idx<=input->rb;idx++)
  {
    if (suftab[idx] > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,suftab[idx]-1,
                                              GT_READMODE_FORWARD);
      if (ISNOTSPECIAL(cc))
      {
        if (gt_bittab_bit_is_set(input->visitor->marktab,cc))
        {
          return false;
        }
        gt_bittab_set_bit(input->visitor->marktab,cc);
      }
    }
  }
  return true;
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
  const GtEncseq *encseq = NULL;
  const GtUword *suftab = NULL;

  gt_assert(ev && err);
  lev = gt_esa_smax_lcpitvs_visitor_cast(ev);
  if (lcp >= lev->searchlength && ret_info->maxlcpinterval)
  {
    lev->input->lb = lb;
    lev->input->rb = rb;
    if (gt_esa_smax_verify_supmax(map_verify_supmax,lev->input))
/*    if (verify_supmax(lev,lb,rb)) */
    {
      for (s=lb;s<rb;s++)
      {
        for (t=s+1;t<=rb;t++)
        {
          if (!lev->silent)
          {
            if (encseq == NULL)
            {
              encseq = gt_encseqSequentialsuffixarrayreader(lev->ssar);
            }
            if (suftab == NULL)
            {
              suftab = gt_suftabSequentialsuffixarrayreader(lev->ssar);
            }
            print_repeat_both(NULL, encseq,
                              lcp, suftab[s], suftab[t],
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
                                      Sequentialsuffixarrayreader *ssar,
                                      GtUword searchlength,
                                      bool absolute,
                                      bool silent,
                                      GT_UNUSED GtProcessSmaxpairs *process_smaxpairs,
                                      GT_UNUSED void *process_smaxpairsdata,
                                      GtTimer *smaxprogress)
{ 
  GtESASmaxLcpintervalsVisitor *lev;
  GtESAVisitor *ev = gt_esa_visitor_create(gt_esa_smax_lcpitvs_visitor_class());
  lev = gt_esa_smax_lcpitvs_visitor_cast(ev);
  lev->ssar = ssar;
  lev->searchlength = searchlength;
  lev->absolute = absolute;
  lev->smaxprogress = smaxprogress;
  lev->silent = silent;
  lev->marktab = gt_bittab_new(gt_alphabet_size(gt_encseq_alphabet(ssar->encseq)));
  lev->input = gt_malloc(sizeof(GtESASmaxMapVerifyInput));
  lev->input->visitor = lev;
  return ev;
}
