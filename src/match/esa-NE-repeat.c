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

Definedunsignedint get_lctx(const GtEncseq *encseq,
                                    GtUword suf)
{
  Definedunsignedint temp;
  if (suf > 0)
  {
/*
    printf("suf: " GT_WU " char: %c\n",suf,
        gt_encseq_get_decoded_char(encseq,
                                   suf,
                                   GT_READMODE_FORWARD));
    printf("suf-1: " GT_WU " char: %c\n",suf-1,
        gt_encseq_get_decoded_char(encseq,
                                   suf-1,
                                   GT_READMODE_FORWARD));
*/
    temp.valueunsignedint = gt_encseq_get_encoded_char(encseq,
                                                suf-1,
                                                GT_READMODE_FORWARD);
    if (ISSPECIAL(temp.valueunsignedint))
    {
      temp.defined = false;
    } else
    {
      temp.defined = true;
    }
  } else {
    temp.defined = false;
  }
  return temp;
}

Definedunsignedint check_lctx(Definedunsignedint lctx_psuf,
                                      Definedunsignedint lctx_nsuf)
{
  if (lctx_psuf.defined) 
  {
    if(lctx_psuf.valueunsignedint == lctx_nsuf.valueunsignedint)
    {    
      lctx_psuf.defined = true;
    } else
    { 
      lctx_psuf.defined = false;
    }
  }
  return lctx_psuf;
}

bool is_notleftextendible(const GtEncseq *encseq,
                              const GtUword *suftab,
                              const GtUword suftab_size)
{
  GtUword s,
          suf,
          suf_comp;
//  GtUword *suftab_copy = gt_malloc(sizeof(GtUword) * suftab_size);
//  seqnum = 0;
/*
  memcpy(suftab_copy, suftab,
         suftab_size * sizeof (GtUword));

  qsort(suftab_copy, suftab_size,
        sizeof (GtUword), &compare_suftabvalues);
*/
  suf = suftab[0];
  if (suf != 0)
  {
    suf_comp = gt_encseq_get_encoded_char(encseq,suf-1,GT_READMODE_FORWARD);
  } else
  {
    suf_comp = 254;
  }
  for (s = 1; s < suftab_size; s++)
  {
    suf = suftab[s];
    if (suf != 0)
    {
/*
      printf(" " GT_WU "",suf);
      printf(" %c",gt_encseq_get_decoded_char(encseq,
                                                     suf-1,
                                                     GT_READMODE_FORWARD));
*/
      if (ISSPECIAL(suf_comp) || ISSPECIAL(gt_encseq_get_encoded_char(encseq,suf-1,GT_READMODE_FORWARD)))
      {
        return true;
      }
      if (suf_comp != gt_encseq_get_encoded_char(encseq,suf-1,GT_READMODE_FORWARD))
      {
        return true;
      }
    } else
    {
      return true;
/*
      printf(" " GT_WU "",suf);
      printf("$");
*/
    }
  }
  return false;
//  printf("\n");
}
