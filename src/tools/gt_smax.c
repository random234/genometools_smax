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
#include "match/esa-scanprj.h"
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
  bool bool_option_sequencedesc;
  bool bool_option_non_extendible;
} SmaxArguments;

typedef struct {
  bool absolute;
  bool direct;
  bool palindromic;
  bool sequencedesc;
  bool compact;
  bool non_extendible;
  char method;
  GtStr *output;
} PrintArguments;

void print_sequence(const GtEncseq *encseq, GtUword pos, GtUword len)
{
  GtUword idx;
  printf("Sequence: ");
  for (idx = pos; idx < pos+len; idx++)
  {
    printf("%c",gt_encseq_get_decoded_char(encseq,idx,GT_READMODE_FORWARD));
  }
  printf("\n");
}

void print_repeat(void *info,
                  const GtEncseq *encseq,
                  GtUword maxlen,
                  GtUword suftab_s,
                  GtUword suftab_t,
                  GtUword seqnum_s,
                  GtUword seqnum_t)
{
  PrintArguments *printinfo = (PrintArguments *) info;
  GtUword score = maxlen * 2;

  if (printinfo->absolute)
  { 
    printf(""GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " " GT_WU"\n",maxlen,
          suftab_s, printinfo->method, maxlen, suftab_t,score);
  } else
  {
    GtUword pos_corr_t = gt_encseq_seqstartpos(encseq, seqnum_t),
                                              pos_corr_s = 
                                              gt_encseq_seqstartpos(
                                              encseq, seqnum_s);
    printf("" GT_WU " " GT_WU " " GT_WU " %3c " GT_WU " " GT_WU " "
           GT_WU " " GT_WU "\n",maxlen, seqnum_s, suftab_s-pos_corr_s,
           printinfo->method, maxlen, seqnum_t, suftab_t-pos_corr_t,score);
  }

  if (printinfo->sequencedesc)
  {
    print_sequence(encseq, seqnum_s, maxlen);
  }
}

void process_repeat_verbose(void *info,
                    const GtEncseq *encseq,
                    const GtUword maxlen,
                    const GtUword *suftab,
                    const GtUword suftab_size)
{
  PrintArguments *printinfo = (PrintArguments *) info;
  GtUword s,t,suftab_s,suftab_t;
  bool reverse_direct;
  for (s = 0;s < suftab_size;s++)
  {
    suftab_s = suftab[s];
    for (t=s+1;t < suftab_size;t++)
    {
      GtUword seqnum_s,seqnum_t;
      suftab_t = suftab[t];
      if (suftab_s > suftab_t)
      {
        GtUword tmp;
        tmp = suftab_s;
        suftab_s = suftab_t;
        suftab_t = tmp;
      }
      seqnum_s = gt_encseq_seqnum(encseq,suftab_s);
      seqnum_t = gt_encseq_seqnum(encseq,suftab_t);
      printinfo->method = 'D';

      if (printinfo->palindromic)
      {
        GtUword halftotal_seqnum = GT_DIV2(gt_encseq_num_of_sequences(encseq));
        if (seqnum_t >= halftotal_seqnum)
        {
          printinfo->method = 'P';
          if (seqnum_s >= halftotal_seqnum)
          {
            /* direct repeat on reverse strand found dropping it */
            reverse_direct = true;
          } 
            suftab_t = GT_REVERSEPOS(gt_encseq_total_length(encseq),suftab_t);
            suftab_t = (suftab_t+1)-maxlen;
            seqnum_t = seqnum_t-halftotal_seqnum;
        }
      }
      if (!reverse_direct)
      {
        print_repeat(info,encseq,maxlen,suftab_s,suftab_t,seqnum_s,seqnum_t);
      } else 
      {
        reverse_direct = false;
      }
      if (printinfo->sequencedesc)
      {
        print_sequence(encseq, suftab_s, maxlen);
      } 
    }
  }
}

void process_repeat_compact(void *info,
                            const GtEncseq *encseq,
                            const GtUword maxlen,
                            const GtUword *suftab,
                            const GtUword suftab_size)
{
  PrintArguments *printinfo = (PrintArguments *) info;
  GtUword s,t,suftab_s,suftab_t;
  bool reverse_direct;

  printf("" GT_WU " (", maxlen);
  for (s = 0;s < suftab_size;s++)
  {
    suftab_s = suftab[s];
    for (t=s+1;t < suftab_size;t++)
    {
      GtUword seqnum_s,seqnum_t;
      suftab_t = suftab[t];
      if (suftab_s > suftab_t)
      {
        GtUword tmp;
        tmp = suftab_s;
        suftab_s = suftab_t;
        suftab_t = tmp;
      }
      seqnum_s = gt_encseq_seqnum(encseq,suftab_s);
      seqnum_t = gt_encseq_seqnum(encseq,suftab_t);
      printinfo->method = 'D';

      if (printinfo->palindromic)
      {
        GtUword halftotal_seqnum = GT_DIV2(gt_encseq_num_of_sequences(encseq));
        if (seqnum_t >= halftotal_seqnum)
        {
          printinfo->method = 'P';
          if (seqnum_s >= halftotal_seqnum)
          {
            /* direct repeat on reverse strand found dropping it */
            reverse_direct = true;
          } 
            suftab_t = GT_REVERSEPOS(gt_encseq_total_length(encseq),suftab_t);
            suftab_t = (suftab_t+1)-maxlen;
            seqnum_t = seqnum_t-halftotal_seqnum;
        }
      }
      if (!reverse_direct)
      {
        if (printinfo->absolute) 
        {
          printf(" " GT_WU " %c", suftab_s, printinfo->method);
          gt_str_append_cstr(printinfo->output," ");
          gt_str_append_ulong(printinfo->output,suftab_t);
          gt_str_append_cstr(printinfo->output," ");
          gt_str_append_char(printinfo->output,printinfo->method);
        } else 
        {
          GtUword pos_corr_t = gt_encseq_seqstartpos(encseq, seqnum_t),
                  pos_corr_s = gt_encseq_seqstartpos(encseq, seqnum_s);
          printf(" " GT_WU " " GT_WU " %c", seqnum_s, suftab_s-pos_corr_s,
                                            printinfo->method);
          gt_str_append_cstr(printinfo->output," ");
          gt_str_append_ulong(printinfo->output,seqnum_t);
          gt_str_append_cstr(printinfo->output," ");
          gt_str_append_ulong(printinfo->output,suftab_t-pos_corr_t);
          gt_str_append_cstr(printinfo->output," ");
          gt_str_append_char(printinfo->output,printinfo->method);
        }
      } else 
      {
        reverse_direct = false;
      }
      if (printinfo->sequencedesc)
      {
        print_sequence(encseq, suftab_s, maxlen);
      } 
    }
  }
  printf(" ) " GT_WU " (%s " GT_WU " )\n", maxlen, gt_str_get(printinfo->output), maxlen*2);
  gt_str_reset(printinfo->output);
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

  option = gt_option_new_bool("s", "Display sequence ",
            &arguments->bool_option_sequencedesc,false);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("ne", "Non-extendible repeats ",
      &arguments->bool_option_non_extendible,false);
  gt_option_parser_add_option(op,option);

  return op;
}

static int gt_smax_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GtError *err)
{
  SmaxArguments *arguments = tool_arguments;
  int had_err = 0;
//  GtFile *file;
//  GtHashmap *prj_hash;
/* write index option parser which checks compatability of supplied indizes */

  gt_error_check(err);
  gt_assert(arguments);

  if (gt_str_array_size(arguments->index_option_smax) == 0) {
    printf("option -ii is mandatory\n");
    had_err = -1;
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
  GtProcessSmaxpairs process_smaxpairs;
  void *process_smaxpairsdata;
  PrintArguments *printargs = gt_malloc(sizeof(PrintArguments));
  int had_err = 0;
  logger = gt_logger_new(arguments->bool_option_smax, GT_LOGGER_DEFLT_PREFIX,
      stdout);

  printargs->absolute = arguments->bool_option_absolute;
  printargs->direct = arguments->bool_option_direct;
  printargs->palindromic = arguments->bool_option_palindromic;
  printargs->sequencedesc = arguments->bool_option_sequencedesc;
  printargs->compact = arguments->bool_option_compact;
  printargs->non_extendible = arguments->bool_option_non_extendible;

  gt_error_check(err);
  gt_assert(arguments);

  /* XXX */
  if (arguments->bool_option_smax)
  {
    printf("argc=%d, parsed_args=%d\n", argc, parsed_args);
  }
  printf("# argv[0]=%s\n", argv[0]);

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
    if (gt_encseq_is_mirrored(gt_encseqSequentialsuffixarrayreader(ssar)) 
        &! printargs->palindromic)
    {
      printf("You supplied a mirrored sequence, please use the -p option for palindromic matches\n");
      had_err = -1;
    } else
    {
      if (printargs->palindromic &!
          gt_encseq_is_mirrored(gt_encseqSequentialsuffixarrayreader(ssar)))
      {
        printf("You requested the -p option, please provide a mirrored sequence\n");
        had_err = -1;
      } else
      {
        if (printargs->compact)
        {
          printargs->output = gt_str_new();
          process_smaxpairsdata = (void *) printargs;
          process_smaxpairs = process_repeat_compact;    
        } else
        {
          process_smaxpairsdata = (void *) printargs;
          process_smaxpairs = process_repeat_verbose;
        }
      }
    }
    if (!had_err)
    {
      if (printargs->non_extendible)
      {
        
      } else 
      {
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
      }
    } 
  } else 
  {
    had_err = -1;
  }

  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  if (printargs->compact)
  {
    gt_str_delete(printargs->output);
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
