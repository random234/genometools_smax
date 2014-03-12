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
#include "match/esa-smaxlcpintervals.h"
#include "match/esa_linsmax.h"
#include "core/error.h"
#include "core/logger.h"

typedef struct {
  bool bool_option_smax;
  GtStrArray  *index_option_smax;
  GtUword ulong_option_searchlength;
  bool bool_option_absolute;
  bool bool_option_silent;
  bool bool_option_linear;
} SmaxArguments;

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
  option = gt_option_new_string_array("esa",
      "Please specify a valid index or indizes",
      arguments->index_option_smax);
  gt_option_parser_add_option(op,option);

  /* GtUword */
  option = gt_option_new_ulong("l", "Please specify the search length",
      &arguments->ulong_option_searchlength,0);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("absolute",
      "Do you want to see absolute position values",
      &arguments->bool_option_absolute,false);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("silent",
      "Start a run without output useful for time measurements",
      &arguments->bool_option_silent,false);
  gt_option_parser_add_option(op,option);

  option = gt_option_new_bool("linear", "Use linear scan implementation",
      &arguments->bool_option_linear,false);
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
    printf("Indexname must not be empty\n");
  }
  return had_err;
}

static int gt_smax_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GtError *err)
{
  SmaxArguments *arguments = tool_arguments;
  GT_UNUSED GtLogger *logger;
  logger = gt_logger_new(arguments->bool_option_smax, GT_LOGGER_DEFLT_PREFIX,
      stdout);

  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* XXX */
  if (arguments->bool_option_smax)
  {
    printf("argc=%d, parsed_args=%d\n", argc, parsed_args);
  }
  printf("# argv[0]=%s\n", argv[0]);

  /*
  gt_runesmaxlcpvalues(const char *inputindex,
                        bool outedges,
                        bool bottomup,
                        GtLogger *logger,
                        GtError *err);
  */
  if (arguments->bool_option_linear)
  {
    gt_runlinsmax(arguments->index_option_smax,
        arguments->ulong_option_searchlength,arguments->bool_option_absolute,
        arguments->bool_option_silent, true, logger, err);

  } else
  {
    gt_runsmaxlcpvalues(arguments->index_option_smax,
        arguments->ulong_option_searchlength,
        arguments->bool_option_absolute,
        arguments->bool_option_silent, true, true, logger, err);
  }
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
