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
#include "core/arraydef.h"
#include "core/bittab_api.h"
#include "match/esa-seqread.h"
#include "match/esa-smax-callbacks.h"

GT_UNUSED static bool gt_smax_verify_supmax(const GtEncseq *encseq,
        const GtArrayGtUlong *suftab_arr, GtBittab *marktab)
{
  GtUword i, array_size;
  GtUword suf;
  gt_bittab_unset(marktab);
        
  array_size = suftab_arr->nextfreeGtUlong;
  for (i = 0;i<array_size;i++)
  { 
    suf = suftab_arr->spaceGtUlong[i];
    if (suf > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,suf-1, 
                                              GT_READMODE_FORWARD);
      if (ISNOTSPECIAL(cc))
      {
        if (gt_bittab_bit_is_set(marktab,cc))
        {
          return false;
        }
        gt_bittab_set_bit(marktab,cc);
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
            pos_corr_s = gt_encseq_seqstartpos( encseq, seqnum_s);
    printf("" GT_WU " " GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " "
        GT_WU " " GT_WU "\n",maxlen, seqnum_s, suftab_s-pos_corr_s, method,
        maxlen, seqnum_t, suftab_t-pos_corr_t,score);
  }
}
