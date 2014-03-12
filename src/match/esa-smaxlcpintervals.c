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

#include "lcpinterval.h"
#include "esa_smax_lcpvisitor.h"
#include "esa-seqread.h"
#include "esa-bottomup.h"
#include "core/encseq_api.h"
#include "core/str_api.h"
#include "esa-dfs.h"
#include "core/showtime.h"

typedef struct  /* global information */
{
  Lcpinterval lastcompletenode;
  int (*processlcpinterval)(void *,const Lcpinterval *);
  void *processinfo;
} Elcpstate;

static Dfsinfo *allocateDfsinfo_elcp(GT_UNUSED Dfsstate *astate)
{
  return (Dfsinfo *) gt_malloc(sizeof (Lcpinterval));
}

static void freeDfsinfo_elcp(Dfsinfo *adfsinfo,GT_UNUSED Dfsstate *state)
{
  gt_free((Lcpinterval *) adfsinfo);
}

static void showbranchingedgeDFS(bool firstsucc,GtUword fd,
                                 GtUword flb,
                                 GtUword sd,GtUword slb)
{
  printf("B %c" GT_WU "" GT_WU "" GT_WU "" GT_WU "\n",firstsucc ? '1' : '0',
      fd,flb,sd,slb);
}

static int processleafedge_elcp(bool firstsucc,
                                GtUword fatherdepth,
                                Dfsinfo *afather,
                                GtUword leafnumber,
                                GT_UNUSED Dfsstate *astate,
                                GT_UNUSED GtError *err)
{
  Lcpinterval *father = (Lcpinterval *) afather;

  printf("L %c" GT_WU "" GT_WU "" GT_WU "\n",firstsucc ? '1' : '0',
                              fatherdepth,father->left,leafnumber);
  return 0;
}

static int processbranchedge_elcp(bool firstsucc,
                                  GtUword fatherdepth,
                                  Dfsinfo *afather,
                                  Dfsinfo *ason,
                                  Dfsstate *astate,
                                  GT_UNUSED GtError *err)
{
  Lcpinterval *father = (Lcpinterval *) afather;
  Lcpinterval *son = (Lcpinterval *) ason;
  Elcpstate *state = (Elcpstate *) astate;

  if (son != NULL)
  {
    showbranchingedgeDFS(firstsucc,fatherdepth,father->left,son->offset,
                       son->left);
  } else
  {
    showbranchingedgeDFS(firstsucc,fatherdepth,father->left,
                       state->lastcompletenode.offset,
                       state->lastcompletenode.left);
  }
  return 0;
}

static int processcompletenode_elcp(
                          GtUword nodeptrdepth,
                          Dfsinfo *anodeptr,
                          GT_UNUSED GtUword nodeptrminusonedepth,
                          Dfsstate *astate,
                          GT_UNUSED GtError *err)
{
  Lcpinterval *nodeptr = (Lcpinterval *) anodeptr;
  Elcpstate *state = (Elcpstate *) astate;

  gt_assert(state != NULL);
  gt_assert(nodeptr != NULL);
  nodeptr->offset = state->lastcompletenode.offset = nodeptrdepth;
  state->lastcompletenode.left = nodeptr->left;
  state->lastcompletenode.right = nodeptr->right;
  if (state->processlcpinterval != NULL)
  {
    if (state->processlcpinterval(state->processinfo,
                                  &state->lastcompletenode) != 0)
    {
      return -1;
    }
  }
  return 0;
}

static void assignleftmostleaf_elcp(Dfsinfo *adfsinfo,
                                    GtUword leftmostleaf,
                                    GT_UNUSED Dfsstate *dfsstate)
{
  ((Lcpinterval *) adfsinfo)->left = leftmostleaf;
}

static void assignrightmostleaf_elcp(Dfsinfo *adfsinfo,
                                     GtUword currentindex,
                                     GT_UNUSED GtUword previoussuffix,
                                     GT_UNUSED GtUword currentlcp,
                                     GT_UNUSED Dfsstate *dfsstate)
{
  ((Lcpinterval *) adfsinfo)->right = currentindex;
}

static int gt_enumlcpvalues(bool outedges,
                            Sequentialsuffixarrayreader *ssar,
                            int (*processlcpinterval)(void *,
                                                      const Lcpinterval *),
                            void *processinfo,
                            GtLogger *logger,
                            GtError *err)
{
  Elcpstate *state;
  bool haserr = false;

  state = gt_malloc(sizeof (*state));
  if (outedges)
  {
    state->processlcpinterval = NULL;
    state->processinfo = NULL;
  } else
  {
    state->processlcpinterval = processlcpinterval;
    state->processinfo = processinfo;
  }
  if (gt_depthfirstesa(ssar,
                       allocateDfsinfo_elcp,
                       freeDfsinfo_elcp,
                       outedges ? processleafedge_elcp : NULL,
                       outedges ? processbranchedge_elcp : NULL,
                       processcompletenode_elcp,
                       assignleftmostleaf_elcp,
                       assignrightmostleaf_elcp,
                       (Dfsstate *) state,
                       logger,
                       err) != 0)
  {
    haserr = true;
  }
  if (!haserr && state->processlcpinterval != NULL)
  {
    state->lastcompletenode.offset = 0;
    state->lastcompletenode.left = 0;
    state->lastcompletenode.right
      = gt_Sequentialsuffixarrayreader_totallength(ssar);
    if (state->processlcpinterval(state->processinfo,
                                  &state->lastcompletenode) != 0)
    {
      haserr = true;
    }
  }
  gt_free(state);
  return haserr ? -1 : 0;
}

static int showlcpinterval(GT_UNUSED void *data,const Lcpinterval *lcpinterval)
{
  printf("N" GT_WU "" GT_WU "" GT_WU "\n",lcpinterval->offset,
                           lcpinterval->left,
                           lcpinterval->right);
  return 0;
}

int gt_runsmaxlcpvalues(GtStrArray *inputindex,
    GtUword searchlength,
    bool absolute,
    bool silent,
    GT_UNUSED bool outedges,
    bool bottomup,
    GtLogger *logger,
    GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;
  GtTimer *smaxprogress = NULL;
  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(gt_str_array_get(
        inputindex,0),
      SARR_LCPTAB |
      SARR_SUFTAB |
      SARR_ESQTAB,
      false,
      /*SEQ_scan,*/
      logger,
      err);

  gt_showtime_enable();
  if (gt_showtime_enabled())
  {
    smaxprogress = gt_timer_new_with_progress_description(
        "finding supermaximal repeats");
    gt_timer_start(smaxprogress);
  }

  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (bottomup)
    {
      GtESAVisitor *elv = gt_esa_smax_lcpitvs_visitor_new(
          ssar, searchlength, absolute, silent, smaxprogress);
      if (gt_esa_bottomup(ssar, elv, err) != 0)
      {
//        printf("gt_esa_bottomup called\n");
        haserr = true;
      } else {
//        printf("gt_esa_bottomup called without errors\n");
      }
      gt_esa_visitor_delete(elv);
    } else
    {
      if (gt_enumlcpvalues(outedges, ssar, showlcpinterval, NULL,
                           logger, err) != 0)
      {
        printf("gt_enumlcpvalues called\n");
        haserr = true;
      }
    }
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }

  if (smaxprogress != NULL)
  {
/* gt_timer_show_progress(smaxprogress,"%GT_WD.%06GT_WDs real %GT_WDs
 * user %GT_WDs system", stdout); */
    gt_timer_show_progress_final(smaxprogress, stdout);
    gt_timer_delete(smaxprogress);
  }

  return haserr ? -1 : 0;
}

static int gt_esa_scantables(Sequentialsuffixarrayreader *ssar,
                             unsigned int mode,GtError *err)
{
  GtUword lcpvalue,
                previoussuffix = 0,
                idx,
                nonspecials,
                sumsuftab = 0,
                sumlcptab = 0;
  bool haserr = false;

  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  if (mode == 1U)
  {
    for (idx = 0; idx < nonspecials; idx++)
    {
      int retval = gt_nextSequentiallcpvalue(&lcpvalue,ssar,err);

      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
      sumlcptab += lcpvalue; /* silly but guarantees that loop is
                                not eliminated by compiler */
      retval = gt_nextSequentialsuftabvalue(&previoussuffix,ssar);
      gt_assert(retval >= 0);
      if (retval == 0)
      {
        break;
      }
      sumsuftab += previoussuffix; /* silly but guarantees that loop is
                                      not eliminated by compiler */
    }
  } else
  {
    if (mode == 2U)
    {
      for (idx = 0; !haserr && idx < nonspecials; idx++)
      {
        SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
        sumlcptab += lcpvalue; /* silly but guarantees that loop is
                                  not eliminated by compiler */
        SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
        sumsuftab += previoussuffix; /* silly but guarantees that loop is
                                        not eliminated by compiler */
      }
    } else
    {
      gt_error_set(err,"illegal mode %u: use 1 or 2",mode);
      haserr = true;
    }
  }
  if (!haserr)
  {
    printf("sumsuftab="GT_WU"\n",sumsuftab);
    printf("sumlcptab="GT_WU"\n",sumlcptab);
  }
  return haserr ? -1 : 0;
}

int gt_smax_runscanesa(const char *inputindex, unsigned int mode,
                  GtLogger *logger,GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(inputindex,
                                                   SARR_LCPTAB |
                                                   SARR_SUFTAB |
                                                   SARR_ESQTAB,
                                                   false,
                                                   logger,
                                                   err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr && gt_esa_scantables(ssar, mode, err) != 0)
  {
    haserr = true;
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
