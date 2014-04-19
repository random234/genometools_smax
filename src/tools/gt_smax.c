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

#include "core/ma.h"
#include "core/unused_api.h"
#include "tools/gt_smax.h"
#include "core/logger.h"
#include "core/error.h"
#include "match/esa-seqread.h"
#include "match/esa-smax.h"
#include "match/esa-smaxlcpintervals.h"
#include "match/esa-smax-scan.h"

typedef struct {
  bool bool_option_smax;
  GtStrArray  *index_option_smax;
  GtUword ulong_option_searchlength;
  bool bool_option_absolute;
  bool bool_option_silent;
  bool bool_option_map;
  bool bool_option_compact;
  bool bool_option_direct;
  bool bool_option_palindromic;
} SmaxArguments;

typedef struct {
  bool absolute;
  bool direct;
  bool palindromic;
} PrintArguments;


void print_repeat(GT_UNUSED void *info,
                  const GtEncseq *encseq,
                  GtUword maxlen,
                  GtUword suftab_s,
                  GtUword suftab_t)
{
  PrintArguments *printargs = (PrintArguments *) info;
  GtUword score = maxlen * 2;
  GtUword seqnum_s = gt_encseq_seqnum(encseq,suftab_s);
  GtUword seqnum_t = gt_encseq_seqnum(encseq,suftab_t);
  char method = 'D';
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
  if (printargs->absolute)
  { 
    printf(""GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " " GT_WU"\n",maxlen,
          suftab_s, method, maxlen, suftab_t,score);
  } else
  {
    GtUword pos_corr_t = gt_encseq_seqstartpos(encseq, seqnum_t),
                                              pos_corr_s = 
                                              gt_encseq_seqstartpos(
                                              encseq, seqnum_s);
    printf("" GT_WU " " GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " "
           GT_WU " " GT_WU "\n",maxlen, seqnum_s, suftab_s-pos_corr_s, method,
           maxlen, seqnum_t, suftab_t-pos_corr_t,score);
  }
}

void print_repeat_mirrored(GT_UNUSED void *info,
                  const GtEncseq *encseq,
                  GtUword maxlen,
                  GtUword suftab_s,
                  GtUword suftab_t)
{
  PrintArguments *printargs = (PrintArguments *) info;
  GtUword score = maxlen * 2;
  GtUword halftotal_seqnum = GT_DIV2(gt_encseq_num_of_sequences(encseq));
  GtUword seqnum_s = gt_encseq_seqnum(encseq,suftab_s);
  GtUword seqnum_t = gt_encseq_seqnum(encseq,suftab_t);
  char method = 'D';

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
  if (seqnum_s < halftotal_seqnum && seqnum_t >= halftotal_seqnum)
  {
    method = 'P';
  /*  printf("" GT_WU "\n",GT_REVERSEPOS(gt_encseq_total_length(encseq), 0));
    printf("" GT_WU "\n", suftab_t);
    printf("" GT_WU "\n", suftab_t-gt_encseq_seqstartpos(encseq, seqnum_t));
    */
    if (printargs->absolute)
    {
      // falsche Berechnung
      suftab_t = (GT_REVERSEPOS(gt_encseq_total_length(encseq), 
                  (suftab_t-gt_encseq_seqstartpos(encseq, seqnum_t))
                  )-maxlen+1)-gt_encseq_seqstartpos(encseq, halftotal_seqnum);
    } else 
    {
      suftab_t = GT_REVERSEPOS(gt_encseq_total_length(encseq),
                (suftab_t-gt_encseq_seqstartpos(encseq, seqnum_t)))-maxlen+1;
    }
  }

  if (printargs->absolute)
  { 
    printf(""GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " " GT_WU"\n",maxlen,
          suftab_s, method, maxlen, suftab_t-gt_encseq_seqstartpos(encseq, halftotal_seqnum),score);
  } else
  {
    GtUword pos_corr_t = gt_encseq_seqstartpos(encseq, seqnum_t);
    GtUword pos_corr_s = gt_encseq_seqstartpos(encseq, seqnum_s);
    printf("" GT_WU " " GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " "
         GT_WU " " GT_WU "\n",maxlen, seqnum_s, suftab_s-pos_corr_s, method,
         maxlen, seqnum_t, suftab_t-pos_corr_t,score);
  }  
 }

static void* gt_smax_arguments_new(void)
{
  SmaxArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->index_option_smax = gt_str_array_new();
  return arguments;
}

static void gt_smax_arguments_delete(void *tool_arguments)
{
  SmaxArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_array_delete(arguments->index_option_smax);
  gt_free(arguments);
}

static GtOptionParser* gt_smax_option_parser_new(void *tool_arguments)
{
  SmaxArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new(
      "[searchlength] [INDEXFILE]",
      "Gt smax searches for supermaximal repeats by detecting local maxima "
      "in a lcp table.");

  /* bool */
  option = gt_option_new_bool("v","be verbose about the output",
      &arguments->bool_option_smax,false);
  gt_option_parser_add_option(op,option);

  /* string */
  option = gt_option_new_string_array("ii",
      "Please specify a valid index or indizes",
      arguments->index_option_smax);
  gt_option_parser_add_option(op,option);

  /* GtUword */
  option = gt_option_new_ulong("l", "Please specify the search length",
      &arguments->ulong_option_searchlength,2);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("absolute",
      "Do you want to see absolute position values",
      &arguments->bool_option_absolute,false);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("silent",
      "Start a run without output useful for time measurements",
      &arguments->bool_option_silent,false);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("map", "Use map scan implementation",
      &arguments->bool_option_map,false);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("compact", "Show compact output ",
      &arguments->bool_option_compact,false);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("d", "Display only matches on forward strand",
      &arguments->bool_option_direct,true);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("p", "Display palindromic matches ",
      &arguments->bool_option_palindromic,false);
  gt_option_parser_add_option(op,option);

  return op;
}

static int gt_smax_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GtError *err)
{
  SmaxArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (gt_str_array_size(arguments->index_option_smax) == 0) {
    printf("option -ii is mandatory\n");
    had_err = 1;
  }
  return had_err;
}

static int gt_smax_runner(int argc,
                          const char **argv,
                          int parsed_args,
                          void *tool_arguments,
                          GtError *err)
{
  SmaxArguments *arguments = tool_arguments;
  GtLogger *logger;
  Sequentialsuffixarrayreader *ssar;
  PrintArguments *printargs = gt_malloc(sizeof(PrintArguments));
  int had_err = 0;
  logger = gt_logger_new(arguments->bool_option_smax, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  printargs->absolute = arguments->bool_option_absolute;
  printargs->direct = arguments->bool_option_direct;
  printargs->palindromic = arguments->bool_option_palindromic;

  gt_error_check(err);
  gt_assert(arguments);

  /* XXX */
  if (arguments->bool_option_smax)
  {
    printf("argc=%d, parsed_args=%d\n", argc, parsed_args);
  }
  printf("# argv[0]=%s\n", argv[0]);

  GtProcessSmaxpairs process_smaxpairs;
  void *process_smaxpairsdata;
    
 
  if((ssar = gt_newSequentialsuffixarrayreaderfromfile(gt_str_array_get(
                                          arguments->index_option_smax,0),
                                          SARR_LCPTAB |
                                          SARR_SUFTAB |
                                          SARR_ESQTAB,
                                          /* scan suftab and lcptab */
                                          /* scan = true */
                                          /* map = false */
                                          !arguments->bool_option_map, 
                                          logger,
                                          err)))
  {
    if (gt_encseq_is_mirrored(gt_encseqSequentialsuffixarrayreader(ssar)))
    {
      process_smaxpairs = print_repeat_mirrored;
      process_smaxpairsdata = (void *) printargs;
    } else
    {
      process_smaxpairs = print_repeat;
      process_smaxpairsdata = (void *) printargs;
    }

    if (arguments->bool_option_map)
    {
      if (gt_runsmaxlcpvalues(ssar,
                              arguments->ulong_option_searchlength,
                              arguments->bool_option_silent, 
                              true, 
                              process_smaxpairs,
                              process_smaxpairsdata,
                              err) != 0)
      {
        had_err = -1;
      }
    } else 
    {
      if (gt_runlinsmax(ssar,
                        arguments->ulong_option_searchlength,
                        arguments->bool_option_silent,
                        process_smaxpairs,
                        process_smaxpairsdata,
                        err) != 0)
      {
        had_err = -1;
      }
    }
  } else 
  {
    had_err = -1;
  }
  gt_free(printargs); 
  gt_logger_delete(logger);
  return had_err;
}

GtTool* gt_smax(void)
{
  return gt_tool_new(gt_smax_arguments_new,
                  gt_smax_arguments_delete,
                  gt_smax_option_parser_new,
                  gt_smax_arguments_check,
                  gt_smax_runner);
}
