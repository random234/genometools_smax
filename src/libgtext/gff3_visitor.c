/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <libgtext/genome_node.h>
#include <libgtext/genome_visitor_rep.h>
#include <libgtext/gff3_output.h>
#include <libgtext/gff3_parser.h>
#include <libgtext/gff3_visitor.h>

struct GFF3Visitor {
  const GenomeVisitor parent_instance;
  bool version_string_shown;
  unsigned long *id_counter;
  Hashtable *genome_feature_to_id_array,
            *genome_feature_to_unique_id_str;
  GenFile *outfp;
};

typedef struct {
  Hashtable *genome_feature_to_id_array;
  const char *id;
} Add_id_info;

typedef struct {
  bool attribute_shown;
  GenFile *outfp;
} ShowAttributeInfo;

#define gff3_visitor_cast(GV)\
        genome_visitor_cast(gff3_visitor_class(), GV)

static void gff3_version_string(GenomeVisitor *gv)
{
  GFF3Visitor *gff3_visitor = gff3_visitor_cast(gv);
  assert(gff3_visitor);
  if (!gff3_visitor->version_string_shown) {
    genfile_xprintf(gff3_visitor->outfp, "%s   %u\n", GFF_VERSION_PREFIX,
                    GFF_VERSION);
    gff3_visitor->version_string_shown = true;
  }
}

static void gff3_visitor_free(GenomeVisitor *gv, Env *env)
{
  GFF3Visitor *gff3_visitor = gff3_visitor_cast(gv);
  assert(gff3_visitor);
  env_ma_free(gff3_visitor->id_counter, env);
  hashtable_delete(gff3_visitor->genome_feature_to_id_array, env);
  hashtable_delete(gff3_visitor->genome_feature_to_unique_id_str, env);
}

static int gff3_visitor_comment(GenomeVisitor *gv, Comment *c, Env *env)
{
  GFF3Visitor *gff3_visitor;
  env_error_check(env);
  gff3_visitor = gff3_visitor_cast(gv);
  assert(gv && c);
  gff3_version_string(gv);
  genfile_xprintf(gff3_visitor->outfp, "#%s\n", comment_get_comment(c));
  return 0;
}

static int add_id(GenomeNode *gn, void *data, Env *env)
{
  Add_id_info *info = (Add_id_info*) data;
  Array *parent_features = NULL;
  env_error_check(env);
  assert(gn && info && info->genome_feature_to_id_array && info->id);
  parent_features = hashtable_get(info->genome_feature_to_id_array, gn);
  if (!parent_features) {
    parent_features = array_new(sizeof (char*), env);
    hashtable_add(info->genome_feature_to_id_array, gn, parent_features, env);
  }
  array_add(parent_features, info->id, env);
  return 0;
}

static int show_attribute(const char *attr_name, const char *attr_value,
                          void *data, Env *env)
{
  ShowAttributeInfo *info = (ShowAttributeInfo*) data;
  env_error_check(env);
  assert(attr_name && attr_value && info);
  if (strcmp(attr_name, ID_STRING) && strcmp(attr_name, PARENT_STRING)) {
    if (info->attribute_shown)
      genfile_xfputc(';', info->outfp);
    else
      info->attribute_shown = true;
    genfile_xprintf(info->outfp, "%s=%s", attr_name, attr_value);
  }
  return 0;
}

static int gff3_show_genome_feature(GenomeNode *gn, void *data, Env *env)
{
  bool part_shown = false;
  GFF3Visitor *gff3_visitor = (GFF3Visitor*) data;
  GenomeFeature *gf = (GenomeFeature*) gn;
  Array *parent_features = NULL;
  ShowAttributeInfo info;
  unsigned long i;
  int had_err = 0;
  Str *id;

  env_error_check(env);
  assert(gn && gf && gff3_visitor);

  /* output leading part */
  gff3_output_leading(gf, gff3_visitor->outfp);

  /* show unique id part of attributes */
  if ((id = hashtable_get(gff3_visitor->genome_feature_to_unique_id_str, gn))) {
    genfile_xprintf(gff3_visitor->outfp, "%s=%s", ID_STRING, str_get(id));
    part_shown = true;
  }

  /* show parent part of attributes */
  parent_features = hashtable_get(gff3_visitor->genome_feature_to_id_array,
                                  gn);
  if (array_size(parent_features)) {
    if (part_shown)
      genfile_xfputc(';', gff3_visitor->outfp);
    genfile_xprintf(gff3_visitor->outfp, "%s=", PARENT_STRING);
    for (i = 0; i < array_size(parent_features); i++) {
      if (i)
        genfile_xfputc(',', gff3_visitor->outfp);
      genfile_xprintf(gff3_visitor->outfp, "%s",
                      *(char**) array_get(parent_features, i));
    }
    part_shown = true;
  }

  /* show missing part of attributes */
  info.attribute_shown = part_shown;
  info.outfp = gff3_visitor->outfp;
  had_err = genome_feature_foreach_attribute(gf, show_attribute, &info, env);

  /* show terminal newline */
  genfile_xfputc('\n', gff3_visitor->outfp);

  return had_err;
}

static int store_ids(GenomeNode *gn, void *data, Env *env)
{
  GFF3Visitor *gff3_visitor = (GFF3Visitor*) data;
  GenomeFeature *gf = (GenomeFeature*) gn;
  GenomeFeatureType type;
  Add_id_info add_id_info;
  int had_err = 0;
  Str *id;

  env_error_check(env);
  assert(gn && gf && gff3_visitor);
  type = genome_feature_get_type(gf);

  if (genome_node_has_children(gn)) {
    /* increase id counter */
    gff3_visitor->id_counter[type]++;

    /* build id string */
    id = str_new_cstr(genome_feature_type_get_cstr(type), env);
    str_append_ulong(id, gff3_visitor->id_counter[type], env);

    /* store (unique) id */
    hashtable_add(gff3_visitor->genome_feature_to_unique_id_str, gn, id, env);

    /* for each child -> store the parent feature in the hash table */
    add_id_info.genome_feature_to_id_array =
      gff3_visitor->genome_feature_to_id_array,
    add_id_info.id = str_get(id);
    had_err = genome_node_traverse_direct_children(gn, &add_id_info, add_id,
                                                   env);
  }
  return had_err;
}

static int gff3_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                       Env *env)
{
  GFF3Visitor *gff3_visitor;
  int had_err;
  env_error_check(env);
  gff3_visitor = gff3_visitor_cast(gv);

  gff3_version_string(gv);

  had_err = genome_node_traverse_children((GenomeNode*) gf, gff3_visitor,
                                          store_ids, true, env);
  if (!had_err) {
    if (genome_node_is_tree((GenomeNode*) gf)) {
      had_err = genome_node_traverse_children((GenomeNode*) gf, gff3_visitor,
                                              gff3_show_genome_feature, true,
                                              env);
    }
    else {
      /* got a DAG -> traverse bin breadth first fashion to make sure that the
         'Parent' attributes are shown in correct order */
      had_err = genome_node_traverse_children_breadth((GenomeNode*) gf,
                                                      gff3_visitor,
                                                      gff3_show_genome_feature,
                                                      true, env);
    }
  }

  /* reset hashtables */
  hashtable_reset(gff3_visitor->genome_feature_to_id_array, env);
  hashtable_reset(gff3_visitor->genome_feature_to_unique_id_str, env);

  /* terminator */
  genfile_xfputs("###\n", gff3_visitor->outfp);

  return had_err;
}

static int gff3_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                        Env *env)
{
  GFF3Visitor *gff3_visitor;
  env_error_check(env);
  gff3_visitor = gff3_visitor_cast(gv);
  assert(gv && sr);
  /* a sequence region has no children */
  assert(!genome_node_has_children((GenomeNode*) sr));

  gff3_version_string(gv);
  genfile_xprintf(gff3_visitor->outfp, "%s   %s %lu %lu\n", GFF_SEQUENCE_REGION,
                  str_get(genome_node_get_seqid((GenomeNode*) sr)),
                  genome_node_get_start((GenomeNode*) sr),
                  genome_node_get_end((GenomeNode*) sr));
  return 0;
}

const GenomeVisitorClass* gff3_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (GFF3Visitor),
                                          gff3_visitor_free,
                                          gff3_visitor_comment,
                                          gff3_visitor_genome_feature,
                                          gff3_visitor_sequence_region,
                                          NULL };
  return &gvc;
}

GenomeVisitor* gff3_visitor_new(GenFile *outfp, Env *env)
{
  GenomeVisitor *gv = genome_visitor_create(gff3_visitor_class(), env);
  GFF3Visitor *gff3_visitor = gff3_visitor_cast(gv);
  gff3_visitor->version_string_shown = false;
  gff3_visitor->id_counter = env_ma_calloc(env,
                                          genome_feature_type_num_of_features(),
                                          sizeof (unsigned long));
  gff3_visitor->genome_feature_to_id_array = hashtable_new(HASH_DIRECT, NULL,
                                                           (FreeFunc)
                                                           array_delete, env);
  gff3_visitor->genome_feature_to_unique_id_str = hashtable_new(HASH_DIRECT,
                                                                NULL,
                                                                (FreeFunc)
                                                                str_delete,
                                                                env);
  gff3_visitor->outfp = outfp;
  return gv;
}
