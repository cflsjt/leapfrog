#!/usr/bin/env python

import sys
import numpy as np
import pysam
import argparse
import collections
import pdb

def parse_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bam')
    parser.add_argument('-a', '--min_coverage', type=int, default=4,
                        help=("Minimal coverage for a family to be "
                              "exported as a peak - the number "
                              "checked is the maximal coverage for "
                              "that set of mapping reads. (default 4)"))
    parser.add_argument('-m', '--min_diff', type=int, default=5,
                        help=("Stringency to determine if a hit maps "
                              "uniquely. This score is the "
                              "minimal allowed difference between the "
                              "scores of the first and second "
                              "hit before a read is assumed to be "
                              "mapped to it's correct, unique, "
                              "location. (default=5)"))
    parser.add_argument('-c', '--trim_cov', type=float, default=3.,
                        help=("Minimal coverage of a read region "
                              "that is considered part of the "
                              "peak. This number is used to trim "
                              "trailing edges from a peak, or, if "
                              "necessary to cut a peak into an "
                              "number of smaller peaks. Note, if you "
                              "set a value below 0, it will be "
                              "interpreted as a fraction of the "
                              "maximal coverage observed for that "
                              "peak. (default: 3)"))
    parser.add_argument('-C', '--min_trim_cov', type=int, default=3,
                        help=("When using a fraction for -c, this "
                              "value sets a lower limit - anything "
                              "below this value will be trimmed "
                              "for sure. (default 3)"))
    parser.add_argument('-q', '--output_nonunique',
                        action='store_true', default=False,
                        help=("also output peaks which are no likely "
                              "to be uniquely mapped"))
    return parser.parse_args(args)


def extract_clusters(sam):
    """
    Extracts all clusters of overlapping reads from a sorted, indexed pysam object.
    Generates a dictionary per cluster.
    Accepts a pysam object.
    Returns a dictionary generator.
    """

    cluster = {"reads": [],
               "reference": "",
               "start": -1,
               "stop": -1}

    for i, read in enumerate(sam.fetch()):
        read_reference = sam.getrname(read.tid)
        read_start = read.pos
        read_stop = read.pos + read.qlen

        # if read overlaps the current cluster
        if (read_reference == cluster["reference"]) and (read_start < cluster["stop"]):

            # add the read to the current cluster
            cluster["reads"].append(read)
            cluster["stop"] = read_stop

        # else read is the start of a new cluster
        else:

            # yield the previous cluster but skip the first blank cluster
            if i > 0:
                yield cluster

            # create a new cluster dictionary based on the current read
            cluster["reads"] = [read]
            cluster["reference"] = read_reference
            cluster["start"] = read_start
            cluster["stop"] = read_stop

    # ensure the final cluster is not skipped
    yield cluster


def extract_references(sam):
    """
    Takes a sam object and returns a cluster-dictionary per reference.
    Accepts a pysam object.
    Returns a dictionary generator.
    """
    for reference in sam.references:
        cluster = {"reads": list(sam.fetch(reference)),
                   "reference": reference}
        yield cluster


def sub_cluster(parent_cluster, read_subset, **kwargs):
    """
    Returns a modified cluster with a subset of the original reads.
    Additional parameters can be added as **kwargs.
    Automatically recalculates start and stop positions based on the subset of reads passed.
    Accepts a dictionary
    Returns a dictionary
    """
    child_cluster = {}

    # avoid passing reference to parent cluster or parent clusters reads
    for key in parent_cluster:
        if key in ("reads", "start", "stop"):
            pass
        else:
            child_cluster[key] = parent_cluster[key]

    # add new attributes passed as kwargs to child cluster
    for key, value in kwargs.items():
        child_cluster[key] = value

    # add the explicitly passed reads to child cluster
    child_cluster["reads"] = read_subset
    child_cluster["start"] = min([read.pos for read in read_subset])
    child_cluster["stop"] = max([(read.pos + read.qlen) for read in read_subset])
    return child_cluster


def split_gaps(cluster_generator):
    """
    Subdivides read-clusters based on gaps between non-overlapping reads.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for parent_cluster in cluster_generator:

        # Dummy cluster for comparison with initial read
        child_cluster = {"start": -1, "stop": -1}

        for i, read in enumerate(parent_cluster["reads"]):
            read_start = read.pos
            read_stop = read.pos + read.qlen

            # if read overlaps the current cluster
            if read_start < child_cluster["stop"]:

                # add the read to the current cluster
                child_cluster["reads"].append(read)
                child_cluster["stop"] = read_stop

            # else read is the start of a new cluster
            else:

                # yield the previous cluster but skip the first dummy cluster
                if i > 0:

                    yield child_cluster

                # create a new cluster dictionary based on the current read
                child_cluster = sub_cluster(parent_cluster, [read])

        # ensure the final cluster is not skipped
        yield child_cluster


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
        orientations = {"forwards": [],
                        "reverse": []}
        for read in parent_cluster["reads"]:
            if read.is_reverse:
                orientations["reverse"].append(read)
            else:
                orientations["forwards"].append(read)

        for orientation, reads in orientations.items():
            if len(reads) > 0:
                child_cluster = sub_cluster(parent_cluster, reads, orientation=orientation)
                yield child_cluster


def filter_unique(cluster_generator, threshold):
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
            if score >= threshold:
                unique_reads.append(read)
            else:
                pass

        if len(unique_reads) > 0:
            child_cluster = sub_cluster(parent_cluster, unique_reads, read_type='UNIQUE')
            yield child_cluster


def filter_depth(cluster_generator, threshold):
    """
    Filters read-clusters based on maximum read depth.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    pass


def map_depth(cluster):
    """
    Appends a numpy array of read depth.
    Accepts a dictionary.
    Returns a dictionary.
    """
    depth = np.zeros((cluster["stop"] - cluster["start"]))
    for read in cluster["reads"]:
        depth[(read.pos - cluster["start"]):(read.pos + read.qlen - cluster["start"])] += 1
    return depth


def group_clusters(cluster_generator, *args):
    """
    Groups cluster-dictionaries by unique combinations of values for an arbitrary number of keys.
    Groups are dictionaries that contain a list of clusters and the key value pairs used to categorise them.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    groups = {}
    for cluster in cluster_generator:
        group = '_'.join([cluster[key] for key in args])
        if group not in groups:
            groups[group] = {"clusters": []}
            for key in args:
                groups[group][key] = cluster[key]
        groups[group]["clusters"].append(cluster)
    for key, values in groups.items():
        yield values


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
    cluster_generator = split_gaps(cluster_generator)
    cluster_generator = split_families(cluster_generator)
    cluster_generator = split_orientation(cluster_generator)
    cluster_generator = split_gaps(cluster_generator)


if __name__ == '__main__':
    main()