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

#include "esa-smax.h"
#include "core/encseq_api.h"

void gt_esa_smax_print(GtESASmaxPrintFunc print_func,
                       void *data,
                       const struct GtEncseq *encseq,
                       GtUword maxlen,
                       GtUword suftab_s,
                       GtUword suftab_t,
                       char method,
                       bool absolute)
{
  if (print_func != NULL)
  {
    print_func(data,encseq,maxlen,suftab_s,suftab_t,method,absolute);
  }
}

bool gt_esa_smax_verify_supmax(GtESASmaxVerifySupmaxFunc verifysupmax_func,
                              void *data)
{
  if (verifysupmax_func != NULL)
  {
    return verifysupmax_func(data);
  }
  /* throw error */
  return false;
}
