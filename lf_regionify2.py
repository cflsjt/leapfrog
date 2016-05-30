#!/usr/bin/env python

import sys
import numpy as np
import pysam
import argparse
import collections
from sklearn.cluster import DBSCAN
# import hdbscan

def parse_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bam')
#    parser.add_argument('-a', '--min_coverage', type=int, default=4,
#                        help=("Minimal coverage for a family to be "
#                              "exported as a peak - the number "
#                              "checked is the maximal coverage for "
#                              "that set of mapping reads. (default 4)"))
    parser.add_argument('-s', '--stringency', type=int, default=5,
                        help=("Stringency to determine if a hit maps "
                              "uniquely. This score is the "
                              "minimal allowed difference between the "
                              "scores of the first and second "
                              "hit before a read is assumed to be "
                              "mapped to it's correct, unique, "
                              "location. (default=5)"))
#    parser.add_argument('-c', '--trim_cov', type=float, default=3.,
#                        help=("Minimal coverage of a read region "
#                              "that is considered part of the "
#                              "peak. This number is used to trim "
#                              "trailing edges from a peak, or, if "
#                              "necessary to cut a peak into an "
#                              "number of smaller peaks. Note, if you "
#                              "set a value below 0, it will be "
#                              "interpreted as a fraction of the "
#                              "maximal coverage observed for that "
#                              "peak. (default: 3)"))
#    parser.add_argument('-C', '--min_trim_cov', type=int, default=3,
#                        help=("When using a fraction for -c, this "
#                              "value sets a lower limit - anything "
#                              "below this value will be trimmed "
#                              "for sure. (default 3)"))
#    parser.add_argument('-q', '--output_nonunique',
#                        action='store_true', default=False,
#                        help=("also output peaks which are no likely "
#                              "to be uniquely mapped"))
    parser.add_argument('--eps', type=int, default=100,
                        help=("When using the DBSCAN method to identify "
                              "read clusters, eps is the minimum distance "
                              "allowable between two points for inclusion "
                              "in the the same neighbourhood"))
    parser.add_argument('--min_tips', type=int, default=5,
                        help=("When using the DBSCAN method to identify "
                              "read clusters, min_tips is the minimum number "
                              "of read tips found in a single neighbourhood "
                              "in order to count as a cluster"))
    return parser.parse_args(args)


# def extract_clusters(sam):
#     """
#     Extracts all clusters of overlapping reads from a sorted, indexed pysam object.
#     Generates a dictionary per cluster.
#     Accepts a pysam object.
#     Returns a dictionary generator.
#     """
#
#     cluster = {"reads": [],
#                "reference": "",
#                "start": -1,
#                "stop": -1}
#
#     for i, read in enumerate(sam.fetch()):
#         read_reference = sam.getrname(read.tid)
#         read_start = read.pos
#         read_stop = read.pos + read.qlen
#
#         # if read overlaps the current cluster
#         if (read_reference == cluster["reference"]) and (read_start < cluster["stop"]):
#
#             # add the read to the current cluster
#             cluster["reads"].append(read)
#             cluster["stop"] = read_stop
#
#         # else read is the start of a new cluster
#         else:
#
#             # yield the previous cluster but skip the first blank cluster
#             if i > 0:
#                 yield cluster
#
#             # create a new cluster dictionary based on the current read
#             cluster["reads"] = [read]
#             cluster["reference"] = read_reference
#             cluster["start"] = read_start
#             cluster["stop"] = read_stop
#
#     # ensure the final cluster is not skipped
#     yield cluster


def extract_references(sam):
    """
    Takes a sam object and returns a cluster-dictionary per reference.
    Accepts a pysam object.
    Returns a dictionary generator.
    """
    for reference, length in zip(sam.references, sam.lengths):
        cluster = {"reads": list(sam.fetch(reference)),
                   "reference": reference,
                   "start": 0,
                   "stop": int(length)}
        yield cluster


def sub_cluster(parent_cluster, read_subset, **kwargs):
    """
    Returns a modified cluster with a subset of the original reads.
    Additional parameters can be added as **kwargs.
    Accepts a dictionary
    Returns a dictionary
    """
    child_cluster = {}

    # avoid passing reference to parent cluster or parent clusters reads
    for key in parent_cluster:
        if key == "reads":
            pass
        else:
            child_cluster[key] = parent_cluster[key]

    # add new attributes passed as kwargs to child cluster
    for key, value in kwargs.items():
        child_cluster[key] = value

    # add the explicitly passed reads to child cluster
    child_cluster["reads"] = read_subset

    return child_cluster


# def split_gaps(cluster_generator):
#     """
#     Subdivides read-clusters based on gaps between non-overlapping reads.
#     Accepts a dictionary generator.
#     Returns a dictionary generator.
#     """
#     for parent_cluster in cluster_generator:
#
#         # Dummy cluster for comparison with initial read
#         child_cluster = {"start": -1, "stop": -1}
#
#         for i, read in enumerate(parent_cluster["reads"]):
#             read_start = read.pos
#             read_stop = read.pos + read.qlen
#
#             # if read overlaps the current cluster
#             if read_start < child_cluster["stop"]:
#
#                 # add the read to the current cluster
#                 child_cluster["reads"].append(read)
#                 child_cluster["stop"] = read_stop
#
#             # else read is the start of a new cluster
#             else:
#
#                 # yield the previous cluster but skip the first dummy cluster
#                 if i > 0:
#
#                     yield child_cluster
#
#                 # create a new cluster dictionary based on the current read
#                 child_cluster = sub_cluster(parent_cluster, [read])
#                 child_cluster["start"] = min([read.pos for read in child_cluster["reads"]])
#                 child_cluster["stop"] = max([(read.pos + read.qlen) for read in child_cluster["reads"]])
#
#         # ensure the final cluster is not skipped
#         yield child_cluster


def split_families(cluster_generator):
    """
    Subdivides read-clusters based on read family.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for parent_cluster in cluster_generator:
        families = collections.defaultdict(list)
        for read in parent_cluster["reads"]:
            try:
                family = read.qname.split('__')[0].rsplit('/', 1)[1]
            except IndexError:
                family = read.qname.split('__')[0]
            families[family].append(read)

        for family, reads in families.items():
            child_cluster = sub_cluster(parent_cluster, reads, family=family)
            yield child_cluster


def split_orientation(cluster_generator):
    """
    Subdivides read-clusters based on read orientation.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for parent_cluster in cluster_generator:
        orientations = {"+": [],
                        "-": []}
        for read in parent_cluster["reads"]:
            if read.is_reverse:
                orientations["-"].append(read)
            else:
                orientations["+"].append(read)

        for orientation, reads in orientations.items():
            if len(reads) > 0:
                child_cluster = sub_cluster(parent_cluster, reads, orientation=orientation)
                yield child_cluster


def filter_unique(cluster_generator, args):
    """
    Filters read-clusters by the amount of uniquely mapped reads.
    threshold=args.min_diff
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for parent_cluster in cluster_generator:
        unique_reads = []
        for read in parent_cluster["reads"]:
            tag_as = -999
            tag_xs = -999
            for tag in read.tags:
                if tag[0] == 'AS':
                    tag_as = tag[1]
                if tag[0] == 'XS':
                    tag_xs = tag[1]

            score = tag_as - tag_xs
            if score >= args.stringency:
                unique_reads.append(read)
            else:
                pass

        if len(unique_reads) > 0:
            child_cluster = sub_cluster(parent_cluster, unique_reads, read_type='UNIQUE')
            yield child_cluster


# def filter_depth(cluster_generator, threshold):
#     """
#     Filters read-clusters based on maximum read depth.
#     Accepts a dictionary generator.
#     Returns a dictionary generator.
#     """
#     pass


def read_depth(cluster):
    """
    Calculates the read depth of a cluster.
    Accepts a dictionary.
    Returns an array.
    """
    depth = np.zeros((cluster["stop"] - cluster["start"]))
    for read in cluster["reads"]:
        depth[(read.pos - cluster["start"]):(read.pos + read.qlen - cluster["start"])] += 1
    return depth


def read_tips(cluster):#
    """
    Returns the read end positions of a cluster based on orientation.
    Returns right tip of forwards reads and left tip of reverse reads
    Accepts a dictionary.
    Returns an array.
    """
    tips = np.zeros(len(cluster["reads"]))
    count = 0
    for read in cluster["reads"]:
        if read.is_reverse:
            tips[count] = read.pos
        else:
            tips[count] = read.pos + read.qlen
        count += 1
    return tips.astype(np.int)


# def group_clusters(cluster_generator, *args):
#     """
#     Groups cluster-dictionaries by unique combinations of values for an arbitrary number of keys.
#     Groups are dictionaries that contain a list of clusters and the key value pairs used to categorise them.
#     Accepts a dictionary generator.
#     Returns a dictionary generator.
#     """
#     groups = {}
#     for cluster in cluster_generator:
#         group = '_'.join([cluster[key] for key in args])
#         if group not in groups:
#             groups[group] = {"clusters": []}
#             for key in args:
#                 groups[group][key] = cluster[key]
#         groups[group]["clusters"].append(cluster)
#     for key, values in groups.items():
#         yield values


# def identify_features_by_std(cluster_generator, *args):
#     """
#     Identifies features by identifying loci with a read depth of two standard deviations above the mean.
#     Clusters are grouped together by a combination attributes as specified by *args.
#     This allows for calculating the mean and standard deviation across multiple references.
#     Accepts a dictionary generator.
#     Returns a dictionary generator.
#     """
#     group_generator = group_clusters(cluster_generator, *args)
#     for group in group_generator:
#         group["depth"] = np.empty(0, dtype = int)
#         for cluster in group["clusters"]:
#             cluster["depth"] = read_depth(cluster)
#             group["depth"] = np.concatenate((group["depth"], cluster["depth"]))
#         group["mean"] = group["depth"].mean()
#         group["std"] = group["depth"].std()
#         group["threshold"] = group["mean"] + (2 * group["std"])
#         for cluster in group["clusters"]:
#             cluster["threshold"] = group["threshold"]
#             cluster["feature"] = cluster["depth"] > cluster["threshold"]
#             yield cluster


# def identify_features_by_cov(cluster_generator, args):
#     """
#     Identifies features that are over a minimum threshold depth in args
#     :param cluster_generator:  a dictionary generator
#     :param args: command line arguments
#     :return: a dictionary generator
#     """
#     threshold = args.trim_cov
#     for cluster in cluster_generator:
#         cluster["depth"] = read_depth(cluster)
#         cluster["feature"] = cluster["depth"] > threshold
#         yield cluster


def identify_features_by_dbscan(cluster_generator, args):
    """
    Identifies features based on a DBSCAN clustering algorithm on tip positions
    :param args: command line arguments
    :param cluster_generator:  a dictionary generator
    :return: a dictionary generator
    """
    for cluster in cluster_generator:
        tips = read_tips(cluster).astype(np.int)
        input_tips = np.array(zip(tips, np.zeros(len(tips))), dtype=np.int)
        dbscan = DBSCAN(eps=args.eps, min_samples=args.min_tips).fit(input_tips)
        cluster["feature"] = np.zeros((cluster["stop"] - cluster["start"]), dtype=bool)
        labels = dbscan.labels_.astype(np.int)
        groups = np.unique(labels)
        groups = groups[groups >= 0]
        for group in groups:
            group_tips = tips[labels == group]
            cluster["feature"][min(group_tips) - 1: max(group_tips)] = True
        cluster["depth"] = read_depth(cluster)
        yield cluster


# def identify_features_by_hdbscan(cluster_generator, min_tips=5):
#     """
#     Identifies features based on a DBSCAN clustering algorithm on tip positions
#     :param eps: maximum distance between read tips to be considered in the same neighbourhood
#     :param min_tips: minimum number of tips in a neighbourhood to count as a feature
#     :param cluster_generator:  a dictionary generator
#     :return: a dictionary generator
#     """
#     for cluster in cluster_generator:
#         if (len(cluster["reads"])) < min_tips:
#             continue
#         tips = read_tips(cluster)
#         input_tips = np.array(zip(tips, np.zeros(len(tips))), dtype=np.int)
#         hdb = hdbscan.HDBSCAN(min_cluster_size=min_tips)
#         cluster["feature"] = np.zeros((cluster["stop"] - cluster["start"]), dtype=bool)
#         labels = hdb.fit_predict(input_tips).astype(np.int)
#         groups = np.unique(labels)
#         groups = groups[groups >= 0]
#         for group in groups:
#             group_tips = tips[labels == group]
#             cluster["feature"][min(group_tips) - 1: max(group_tips)] = True
#         cluster["depth"] = read_depth(cluster)
#         yield cluster


def extract_features(cluster_generator):
    """
    Extracts features for a gff file from clusters based of the feature attribute-array.
    Features are returned as dictionaries with start and stop attributes and inherit other attributes
    from their parent cluster.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for cluster in cluster_generator:

        feature = {"start": 0, "stop": 0, "mean_depth": 0}
        for key in cluster:
            if key not in ["reads", "start", "stop", "feature", "depth"]:
                feature[key] = cluster[key]
            else:
                pass

        # determine position of features
        previously_in_feature = False
        for position, currently_in_feature in enumerate(cluster["feature"]):

            if not previously_in_feature and currently_in_feature:
                # start of a feature
                previously_in_feature = True
                feature["start"] = position + 1

            elif previously_in_feature and not currently_in_feature:
                # end of a feature
                previously_in_feature = False
                feature["stop"] = position
                feature["mean_depth"] = cluster["depth"][feature["start"] - 1: feature["stop"]].mean()
                feature["depth"] = cluster["depth"][feature["start"] - 1: feature["stop"]]
                yield feature


def format_features(feature_generator):
    """
    Formats feature dictionaries as strings for a gff file.
    Accepts a dictionary generator.
    Returns a string generator.
    """
    for feature in feature_generator:
        formated = "\t".join([feature["reference"],
                              "REFS",
                              "REFS." + feature["read_type"] + "." + feature["family"],
                              str(feature["start"]),
                              str(feature["stop"]),
                              str(feature["mean_depth"]),
                              feature["orientation"],
                              ".",
                              "ID=reps" + feature["reference"] + str(feature["start"]) +
                              feature["family"] + ";Name=" + feature["family"]])
        yield formated


def output_features(formatted_features):
    for feature in formatted_features:
        print(feature)


def trim_clusters():
    """
    Trims read-clusters based on a read depth.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    pass


def construct_gff():
    """
    Accepts a dictionary generator.
    Prints a gff file to standard out.
    """
    pass


def main():
    args = parse_args(sys.argv[1:])
    sam = pysam.Samfile(args.input_bam, 'rb')
    cluster_generator = extract_references(sam)
    cluster_generator = split_families(cluster_generator)
    cluster_generator = split_orientation(cluster_generator)
    cluster_generator = filter_unique(cluster_generator, args)
    cluster_generator = identify_features_by_dbscan(cluster_generator, args)
    feature_generator = extract_features(cluster_generator)
    formatted_features = format_features(feature_generator)
    output_features(formatted_features)


if __name__ == '__main__':
    main()