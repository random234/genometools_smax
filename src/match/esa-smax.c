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

#include "core/unused_api.h"
#include "core/encseq_api.h"
#include "core/arraydef.h"
#include "core/bittab_api.h"
#include "match/esa-seqread.h"
#include "match/esa-smax.h"

bool gt_esa_smax_verify_supmax(const GtEncseq *encseq,
                              const GtUword *suftabsubset,
                              const GtUword suftabsubset_size,
                              GtBittab *marktab)
{
  GtUword idx;
  gt_bittab_unset(marktab);

  for (idx=0;idx<suftabsubset_size;idx++)
  {
    GtUword suf = suftabsubset[idx];
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
