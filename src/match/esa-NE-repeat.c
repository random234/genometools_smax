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
#include "core/stack-inlined.h"
#include "match/esa-seqread.h"
#include "match/esa-smax.h"
#include "core/defined-types.h"

/* vorraussichtlice Stackgroesse muss geschaetzt werden um
 * konstante Stackgroesse sinnvoll zu waehlen.
 */

typedef struct
{
  GtUword lcp;
  GtUword lb;
  Definedunsignedint left_context;
  //char left_context;
  //bool left_extendible;
} Lcp_stackelem;

GT_STACK_DECLARESTRUCT(Lcp_stackelem, 32UL);


/*
char LEletter(char bwt1, char bwt2)
{
  if (bwt1 == '$' || bwt1 != bwt2)
  {
    return '$';
  } else
  {
    return bwt1;
  }
}
*/

bool check_suffix_not_null(GtUword suffix)
{
  if (suffix == 0) 
  {
    return false;
  } else
  {
    return true;
  }
}

Definedunsignedint check_suf_extendible(const GtEncseq *encseq, 
                                          Definedunsignedint suf1, 
                                          Definedunsignedint suf2)
{
  if (!suf1.defined) 
  {
    return suf1;
  }
  char sufchar1 = gt_encseq_get_decoded_char(encseq,
                                             suf1.valueunsignedint,
                                             GT_READMODE_FORWARD);
  char sufchar2 = gt_encseq_get_decoded_char(encseq,
                                              suf2.valueunsignedint,
                                              GT_READMODE_FORWARD);
  if (ISSPECIAL(sufchar1) || ISSPECIAL(sufchar2) || (sufchar1 != sufchar2))
  {
    suf1.defined = false;
    return suf1;
  } else 
  {
    return suf1;
  }
}

int gt_run_NE_repeats(Sequentialsuffixarrayreader *ssar,
                      GT_UNUSED GtUword searchlength,
                      GT_UNUSED bool silent,
                      GT_UNUSED GtProcessSmaxpairs process_smaxpairs,
                      GT_UNUSED void *process_smaxpairsdata,
                      GtError *err)
{
  bool haserr = false;
  gt_error_check(err);
  GtTimer *nerepeat_progress = NULL;
  GtStackLcp_stackelem lcpstack;
  Lcp_stackelem current_elem;
  GT_UNUSED GtUword lcpvalue = 0,
          previoussuffix = 0,
          idx,
          nonspecials,
          lb;
  Definedunsignedint suf,
          suf_prev,
          suf_curr;
          

  const GtEncseq *encseq = gt_encseqSequentialsuffixarrayreader(ssar);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  GT_STACK_INIT(&lcpstack,32UL);

/*  gt_showtime_enable();
  if (gt_showtime_enabled())
  {
    nerepeat_progress = gt_timer_new_with_progress_description(
        "finding non-extendible repeats");
  gt_timer_start(nerepeat_progress);
  }
*/

//  bwt1 = gt_encseq_get_encoded_char(encseq, 0, GT_READMODE_FORWARD);

/*  
  current_elem.lcp = 0;
  current_elem.lb = 0;
  current_elem.left_context = suf_pred;
  current_elem.left_extendible = false;
*/

//  GT_STACK_PUSH(&lcpstack,current_elem);
  for (idx = 0; idx < nonspecials; idx++)
  {
    if (idx == 0)
    {
      lb = idx;
      SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
      current_elem.lcp = 0;
      current_elem.lb = lb;
      if (check_suffix_not_null(previoussuffix))
      {
        suf_prev.valueunsignedint = previoussuffix-1;
        suf_prev.defined = true;
        current_elem.left_context = suf_prev;
      } else
      {
        suf_prev.defined = false;
        current_elem.left_context = suf_prev;
      }
      GT_STACK_PUSH(&lcpstack,current_elem);
    } else
    {
      lb = idx;
      SSAR_NEXTSEQUENTIALLCPTABVALUE(lcpvalue,ssar);
      SSAR_NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);

      if (suf_prev.defined)
      {
        if (check_suffix_not_null(previoussuffix))
        {
          suf_curr.valueunsignedint = previoussuffix-1;
          suf_curr.defined = true;
        }
        suf = check_suf_extendible(encseq, suf_prev, suf_curr);
      } else
      {
        suf = suf_prev;
      }

      while (GT_STACK_TOP(&lcpstack).lcp > lcpvalue)
      {
        current_elem = GT_STACK_POP(&lcpstack);
        if (!current_elem.left_context.defined && lcpvalue > searchlength)
        {
          if (!silent)
          {
            printf("" GT_WU " " GT_WU " " GT_WU "\n",current_elem.lcp,
                                                current_elem.lb, lb);
          }
        }
        lb = current_elem.lb;
        GT_STACK_TOP(&lcpstack).left_context = check_suf_extendible(
                                      encseq,
                                      current_elem.left_context,
                                      GT_STACK_TOP(&lcpstack).left_context);
        suf = check_suf_extendible(encseq, current_elem.left_context,suf);
      }

      if (GT_STACK_TOP(&lcpstack).lcp == lcpvalue)
      {
        GT_STACK_TOP(&lcpstack).left_context = check_suf_extendible(
                                      encseq,
                                      GT_STACK_TOP(&lcpstack).left_context,
                                      suf);
      } else
      {
        current_elem.lcp = lcpvalue;
        current_elem.lb = lb;
        current_elem.left_context = suf;
        GT_STACK_PUSH(&lcpstack,current_elem);
      }
    }
//      printf("lcp: " GT_WU "\t",lcpvalue);
//      printf("suf: " GT_WU "\t",previoussuffix);
//      if (check_suffix_extendible(previoussuffix))
//      {
//        printf("char: %c\n",  gt_encseq_get_decoded_char(encseq, previoussuffix-1,GT_READMODE_FORWARD));
//      }
//    }
//    suf_pred = gt_encseq_get_decoded_char(encseq,previoussuffix-1,GT_READMODE_FORWARD);
//    suf_curr = gt_encseq_get_decoded_char(encseq,previoussuffix, GT_READMODE_FORWARD);
//    suf = 0; 
/*
    while (GT_STACK_TOP(&lcpstack).lcp > lcpvalue)
    {
      current_elem = GT_STACK_POP(&lcpstack);
      if (current_elem.bwt == '$' && lcpvalue > searchlength)
      {
        if (!silent)
        {
          
          process_smaxpairs (process_smaxpairsdata,
                              encseq,
                              current_elem.lcp,
                              suftab_arr.spaceGtUword,
                              suftab_arr.nextfreeGtUword);
          
          printf("" GT_WU " " GT_WU " " GT_WU "\n",
                current_elem.lcp,current_elem.lb,idx);
        }
      }
      lb = current_elem.lb;
      GT_STACK_TOP(&lcpstack).bwt = LEletter(current_elem.bwt,
                                            GT_STACK_TOP(&lcpstack).bwt);
      bwt = LEletter(current_elem.bwt,bwt);
    }
    if (GT_STACK_TOP(&lcpstack).lcp == lcpvalue)
    {
      GT_STACK_TOP(&lcpstack).bwt = LEletter(GT_STACK_TOP(&lcpstack).bwt, bwt);
    } else
    {
      current_elem.lcp = lcpvalue;
      current_elem.lb = lb;
      current_elem.bwt = bwt;
      GT_STACK_PUSH(&lcpstack, current_elem);
    }
    */
  }

  if (nerepeat_progress != NULL)
  {
    gt_timer_show_progress_final(nerepeat_progress, stdout);
    gt_timer_delete(nerepeat_progress);
  }

/*  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
*/
  return haserr ? -1 : 0;
}
