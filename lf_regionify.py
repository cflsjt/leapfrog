#!/usr/bin/env python

import sys
import numpy as np
import pysam
import argparse
import collections


def count_unique_reads(reads, args):
    unique_reads = 0
    for read in reads:
        AS = -999
        XS = -999
        SDF = 0
        for t in read.tags:
            if t[0] == 'AS':
                AS = t[1]
            if t[0] == 'XS':
                XS = t[1]

        SDF = AS - XS
        # lg.critical("%s %s %s" % ( AS, XS, SDF))
        if SDF >= args.min_diff:
            unique_reads += 1
    return unique_reads


def gap_split(reads, cluster_start, cluster_map, chromosome, family, args):
    gapsplit_rv = 0

    # print "SPLIT"
    # print 'start, stop, maxcov', start, stop, maxcov
    # print covarr
    # after trimming start & stop we should always start inside a
    # peak
    def _process_peak(peakstart, peakstop):
        """
        process a peak - (iteratively send it of for processing

        """
        # identify reads that overlap with this peak
        peakreads = []
        for r in reads:
            # print r.pos, r.pos + r.qlen, start + peakstart, stop + peakstop
            if r.pos > cluster_start + peakstop:
                continue
            if r.pos + r.qlen < cluster_start + peakstart:
                continue
            peakreads.append(r)
        # print 'peak %d to %d with %d reads out of %d ' % (
        # peakstart, peakstop, len(peakreads), len(reads))
        return famdump(chromosome, family, peakreads, args)

    in_peak, peak_start, peak_stop = True, 0, 0
    for current_position, below_threshold in enumerate(cluster_map == 0):
        if not below_threshold:
            # we're not in a gap (yet)
            peak_stop = current_position
            if not in_peak:
                # and we're just coming into a peak
                in_peak = True
                peak_start = current_position
        else:
            # we are in a gap -
            if in_peak:
                # we just left a peak
                gapsplit_rv += _process_peak(peak_start, peak_stop)
                cluster_map[peak_start:peak_stop+1] = 99
            in_peak = False
    # we're leaving a peak for sure now
    if in_peak:
        gapsplit_rv += _process_peak(peak_start, peak_stop)
        cluster_map[peak_start:peak_stop+1] = 99

    return gapsplit_rv


def famdump(chromosome, family, reads, args):

    unique_reads = count_unique_reads(reads, args)

    # calculate the bounds of the clustered reads
    cluster_start = min([read.pos for read in reads])
    cluster_stop = max([read.pos + read.qlen for read in reads])

    # generate a depth map of the clustered reads
    cluster_map = np.zeros((cluster_stop - cluster_start))
    for read in reads:
        cluster_map[(read.pos - cluster_start):(read.pos + read.qlen - cluster_start)] += 1

    # max depth of clustered reads
    cluster_depth = max(cluster_map)

    # set any read depths less than args.trim_cov to 0
    # thresholds less than 1 are treated as a proportion
    threshold_depth = args.trim_cov
    if threshold_depth < 1:
        threshold_depth *= cluster_depth
        threshold_depth = max(args.min_trim_cov, threshold_depth)

    cluster_map[cluster_map < threshold_depth] = 0

    # trim read depths of 0 (less than the threshold) from the ends of the cluster map
    # depths of 0 between taller peaks must be kept
    while len(cluster_map) > 0 and cluster_map[0] == 0:
        cluster_start += 1
        cluster_map = cluster_map[1:]
    while len(cluster_map) > 0 and cluster_map[-1] == 0:
        cluster_stop -= 1
        cluster_map = cluster_map[:-1]

    # Do not include cluster if no points are above the minimum depth threshold
    if len(cluster_map) == 0:
        return 0
    if cluster_depth < args.min_coverage:
        return 0

    # sum of length of individual reads / length of cluster
    mean_read_depth = sum([read.qlen for read in reads]) / float(cluster_stop - cluster_start)

    if 0 in cluster_map:
        # go into gap splitting mode
        return gap_split(reads, cluster_start, cluster_map, chromosome, family, args)

    # Munging read names
    names = [read.qname.split('__')[0].rsplit('/', 1)[0] for read in reads]
    name_count = reversed(sorted([(float(names.count(e)) / len(names), e) for e in set(names)]))
    name_count = [x for x in name_count if x[0] > 0.1]
    names = ",".join([n[1] for n in name_count])

    # Identify direction of reads in cluster
    proportion_reversed = float(len([1 for read in reads if read.is_reverse])) / len(reads)
    if proportion_reversed < 0.1:
        strand = '+'
    elif proportion_reversed < 0.9:
        strand = '.'
    else:
        strand = '-'

    # Only output read cluster if it contains enough unique reads
    if unique_reads >= args.min_diff:
        # Carry on because read cluster contains enough unique reads
        cluster_type = 'UNIQUE'
    elif args.output_nonunique:
        # Carry on regardless of unique read count
        cluster_type = 'NOTUNIQUE'
    else:
        # Not enough unique reads in read cluster
        return 0

    # Print the valid read cluster
    print('%s\tREFS\tREFS.%s.%s\t%d\t%d\t%.2f\t%s\t.\tID=reps_%s_%s_%s;Name="%s"' % (
        chromosome, cluster_type, family, cluster_start, cluster_stop, mean_read_depth,
        strand, chromosome, cluster_start, family, names))

    return 1


def bumpdump(chrom, reads, args):

    fams = collections.defaultdict(list)

    # first subdivide into families
    for read in reads:
        try:
            family = read.qname.split('__')[0].rsplit('/', 1)[1]
        except IndexError:
            family = read.qname.split('__')[0]

        fams[family].append(read)

    rv = 0
    for f in fams:
        rv += famdump(chrom, f, fams[f], args)
    return rv


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


def regionify(sam, args):
    bump_chrom = ""
    bump_start = -1
    bump_stop = -1
    bump_reads = []
    bump_count = 0
    for i, read in enumerate(sam.fetch()):
        chrom = sam.getrname(read.tid)
        start = read.pos
        stop = read.pos + read.qlen
        if bump_chrom != chrom or start > bump_stop:
            # new bump - process old bump
            if bump_reads:
                bump_count += bumpdump(bump_chrom, bump_reads, args)
            bump_chrom = chrom
            bump_start = start
            bump_stop = stop
            bump_reads = [read]
        else:
            # still in bump
            # we demand that the reads are sorted - so only need to check the
            bump_stop = max(stop, bump_stop)
            bump_reads.append(read)
    bump_count += bumpdump(bump_chrom, bump_reads, args)


def main():
    args = parse_args(sys.argv[1:])
    sam = pysam.Samfile(args.input_bam, 'rb')
    regionify(sam, args)


if __name__ == '__main__':
    main()
