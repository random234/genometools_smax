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
  GtUword suf;
  GtUword lb;
  Definedunsignedint left_context;
} Lcp_stackelem;

GT_STACK_DECLARESTRUCT(Lcp_stackelem, 32UL);

Definedunsignedint get_left_context(const GtEncseq *encseq,
                                    GtUword suf)
{
  Definedunsignedint temp;
  if (suf > 0)
  {
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

Definedunsignedint check_left_context(Definedunsignedint bwt_psuf,
                                      Definedunsignedint bwt_nsuf)
{
  if (bwt_psuf.defined) 
  {
    if(bwt_psuf.valueunsignedint == bwt_nsuf.valueunsignedint)
    {    
      bwt_psuf.defined = true;
    } else
    { 
      bwt_psuf.defined = false;
    }
  }
  return bwt_psuf;
}

