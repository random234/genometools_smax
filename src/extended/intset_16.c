/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

/*
  THIS FILE IS GENERATED by
  ./scripts/gen-intsets.rb.
  DO NOT EDIT.
*/

#include <inttypes.h>
#include <limits.h>

#include "core/assert_api.h"
#include "core/intbits.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "extended/intset_rep.h"
#include "extended/intset_16.h"

#define gt_intset_16_cast(cvar) \
        gt_intset_cast(gt_intset_16_class(), cvar)

struct GtIntset16 {
  GtIntset parent_instance;
  uint16_t *elements;
};

static bool gt_intset_16_binarysearch_is_member(const uint16_t *leftptr,
                                                const uint16_t *rightptr,
                                                uint16_t elem)
{
  const uint16_t *midptr;
    while (leftptr <= rightptr) {
      midptr = leftptr + (((GtUword) (rightptr-leftptr)) >> 1);
      if (elem < *midptr) {
        rightptr = midptr-1;
      }
      else {
        if (elem > *midptr)
          leftptr = midptr + 1;
        else
          return true;
      }
    }
  return false;
}

static GtUword gt_intset_16_binarysearch_pos2seqnum(const uint16_t *leftptr,
                                                    const uint16_t *rightptr,
                                                    uint16_t pos)
{
  const uint16_t *midptr,
                   *found = NULL,
                   *leftorig = leftptr;

  gt_assert(leftptr <= rightptr);
  if (pos < *leftptr)
    return 0;

  if (pos > *rightptr)
    return 1UL + (GtUword) (rightptr - leftptr);

  gt_assert(pos > *leftptr && pos < *rightptr);
  while (leftptr <= rightptr) {
    midptr = leftptr + ((GtUword) (rightptr - leftptr) >> 1);
    if (pos < *midptr)
      rightptr = midptr - 1;
    else {
      found = midptr;
      if (pos > *midptr)
        leftptr = midptr + 1;
      else
        break;
    }
  }
  gt_assert(found != NULL && found >= leftorig);
  return 1UL + (GtUword) (found - leftorig);
}

void gt_intset_16_add(GtIntset *intset, GtUword elem)
{
  GtIntset16 *intset_16 = gt_intset_16_cast(intset);
  GtIntsetMembers *members = intset->members;
  GtUword *secstart = members->sectionstart;
  gt_assert(members->nextfree < members->num_of_elems &&
            elem <= members->maxelement &&
            (members->previouselem == ULONG_MAX ||
                                      members->previouselem < elem));
  while (elem >= GT_SECTIONMINELEM(members->currentsectionnum + 1)) {
    gt_assert(members->currentsectionnum < members->numofsections);
    secstart[members->currentsectionnum + 1] = members->nextfree;
    members->currentsectionnum++;
  }
  gt_assert(GT_SECTIONMINELEM(members->currentsectionnum) <= elem &&
            elem < GT_SECTIONMINELEM(members->currentsectionnum+1) &&
            GT_ELEM2SECTION(members->logsectionsize,elem) ==
            members->currentsectionnum);
  intset_16->elements[members->nextfree++] = (uint16_t) elem;
  members->previouselem = elem;
}

bool gt_intset_16_is_member(GtIntset *intset, GtUword elem)
{
  GtIntset16 *intset_16 = gt_intset_16_cast(intset);
  GtIntsetMembers *members = intset->members;
  GtUword *secstart = members->sectionstart;
  if (elem <= members->maxelement)
  {
    const GtUword sectionnum = GT_ELEM2SECTION(members->logsectionsize, elem);

    if (secstart[sectionnum] < secstart[sectionnum+1]) {
      return gt_intset_16_binarysearch_is_member(
                              intset_16->elements + secstart[sectionnum],
                              intset_16->elements + secstart[sectionnum+1] - 1,
                              (uint64_t) elem);
    }
  }
  return false;
}

GtUword gt_intset_16_pos2seqnum(GtIntset *intset, GtUword pos)
{
  GtIntset16 *intset_16 = gt_intset_16_cast(intset);
  GtIntsetMembers *members = intset->members;

  GtUword sectionnum = GT_ELEM2SECTION(members->logsectionsize,pos);

  gt_assert(pos <= members->maxelement);
  if (members->sectionstart[sectionnum] < members->sectionstart[sectionnum+1]) {
    return members->sectionstart[sectionnum] +
           gt_intset_16_binarysearch_pos2seqnum(
                  intset_16->elements + members->sectionstart[sectionnum],
                  intset_16->elements + members->sectionstart[sectionnum+1] - 1,
                  (uint16_t) pos);
  }
  return members->sectionstart[sectionnum];
}

size_t gt_intset_16_size(GtUword maxelement, GtUword num_of_elems)
{
  size_t logsectionsize = (sizeof (uint16_t)) + CHAR_BIT;
  return sizeof (uint16_t) * num_of_elems +
    sizeof (GtUword) * (GT_ELEM2SECTION(logsectionsize, maxelement) + 1);
}

void gt_intset_16_delete(GtIntset *intset)
{
  GtIntset16 *intset_16 = gt_intset_16_cast(intset);
  if (intset_16 != NULL) {
    gt_free(intset_16->elements);
  }
}

/* map static local methods to interface */
const GtIntsetClass* gt_intset_16_class(void)
{
  static const GtIntsetClass *this_c = NULL;
  if (this_c == NULL) {
    this_c = gt_intset_class_new(sizeof (GtIntset16),
                                 gt_intset_16_add,
                                 gt_intset_16_is_member,
                                 gt_intset_16_pos2seqnum,
                                 gt_intset_16_delete);
  }
  return this_c;
}

GtIntset* gt_intset_16_new(GtUword maxelement, GtUword num_of_elems)
{
  GtIntset *intset;
  GtIntset16 *intset_16;
  GtIntsetMembers *members;
  GtUword idx;

  intset = gt_intset_create(gt_intset_16_class());
  intset_16 = gt_intset_16_cast(intset);
  members = intset->members;

  intset_16->elements =
    gt_malloc(sizeof (*intset_16->elements) * num_of_elems);
  members->logsectionsize = sizeof (uint16_t) * CHAR_BIT;
  members->nextfree = 0;
  members->numofsections = GT_ELEM2SECTION(members->logsectionsize,
                                           maxelement) + 1;
  members->sectionstart = gt_malloc(sizeof (*members->sectionstart) *
                                    (members->numofsections + 1));
  members->sectionstart[0] = 0;
  for (idx = (GtUword) 1; idx <= members->numofsections; idx++) {
    members->sectionstart[idx] = num_of_elems;
  }
  members->maxelement = maxelement;
  members->currentsectionnum = 0;
  members->num_of_elems = num_of_elems;
  members->previouselem = ULONG_MAX;
  return intset;
}

int gt_intset_16_unit_test(GtError *err) {
  int had_err = 0;
  GtIntset *is;
  GtUword num_of_elems = gt_rand_max(((GtUword) 1) << 10) + 1,
          *arr = gt_malloc(sizeof (*arr) * num_of_elems),
          stepsize = (num_of_elems <<4 / num_of_elems) >> 1,
          idx;
  size_t is_size;

  gt_error_check(err);

  arr[0] = gt_rand_max(stepsize) + 1;
  for (idx = (GtUword) 1; idx < num_of_elems; ++idx) {
    arr[idx] = arr[idx - 1] + gt_rand_max(stepsize) + 1;
  }

  is_size = gt_intset_16_size(arr[num_of_elems - 1], num_of_elems);

  if (!had_err) {
    if (is_size < (size_t) UINT_MAX) {
      is = gt_intset_16_new(arr[num_of_elems - 1], num_of_elems);
      for (idx = 0; idx < num_of_elems; idx++) {
        gt_intset_16_add(is, arr[idx]);
      }
      had_err = gt_intset_unit_test_notinset(is, 0, arr[0] - 1, err);
      if (!had_err)
        had_err = gt_intset_unit_test_check_seqnum(is, 0, arr[0] - 1, 0, err);
      for (idx = (GtUword) 1; !had_err && idx < num_of_elems; idx++) {
        had_err = gt_intset_unit_test_notinset(is, arr[idx - 1] + 1,
                                               arr[idx] - 1, err);
        if (!had_err)
          had_err = gt_intset_unit_test_check_seqnum(is, arr[idx - 1] + 1,
                                                     arr[idx] - 1, idx, err);
      }
      gt_intset_16_delete(is);
    }
  }
  gt_free(arr);
  return had_err;
}
