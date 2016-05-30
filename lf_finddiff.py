#!/usr/bin/env python

import os
import string
import argparse
import collections


class Peekorator(object):

    def __init__(self, generator):
        self.empty = False
        self.peek = None
        self.generator = generator
        try:
            self.peek = self.generator.next()
        except StopIteration:
            self.empty = True

    def __iter__(self):
        return self

    def next(self):
        """
        Return the self.peek element, or raise StopIteration
        if empty
        """
        if self.empty:
            raise StopIteration()
        to_return = self.peek
        try:
            self.peek = self.generator.next()
        except StopIteration:
            self.peek = None
            self.empty = True
        return to_return


class GFFRecord:
    def __init__(self, line, tag):
        ls = line.split()
        self.tag = tag
        self.seqid = ls[0]
        self.source = ls[1]
        self.type = ls[2]
        self.start = int(ls[3])
        self.end = int(ls[4])
        self.score = float(ls[5])
        self.strand = ls[6]
        self.phase = ls[7]
        attrs = ls[8].split(';')
        for a in attrs:
            k, v = map(string.strip, a.split('=', 1))
            self.__dict__[k] = v

    def __str__(self):
        return self.ID


class GFFReader(object):
    def __init__(self, filename, only_unique):
        self.filename = filename
        self.F = open(self.filename)
        self.tag = filename.rsplit('/', 1)[-1].replace('.gff', '')
        self.only_unique = only_unique

    def __iter__(self):
        return self

    def getline(self):
        while True:
            line = self.F.readline()
            if not line:
                break
            if line[0] != '#':
                break
        if not line:
            raise StopIteration()
        return line

    def next(self):
        while True:
            record = GFFRecord(self.getline(), self.tag)
            if self.only_unique and record.type[:11] != 'REFS.UNIQUE':
                continue
            break
        return record

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--base')
parser.add_argument('inputgff', nargs='+')
parser.add_argument('-q', '--only_unique', action='store_true', default=False,
                    help='only process "REFS.UNQIUE.*" features')
parser.add_argument('-d', '--only_differential', action='store_true', default=False,
                    help='show only regions that are differential between the samples')

args = parser.parse_args()

base = args.base
inputfiles = args.inputgff
nicenames = [os.path.basename(x).replace('.gff', '') for x in inputfiles]
# outgff = ['%s.%s.gff' % (base, x) for x in nicenames]
parsers = [Peekorator(GFFReader(x, args.only_unique)) for x in inputfiles]

FOUT1 = open(base + '.regions', 'w')
FOUT2 = open(base + '.table', 'w')

COREGFF = open(base + '.gff', 'w')
# GFFOUT = [open(x, 'w') for x in outgff]

# write FOUT2 header
FOUT2.write("\t")
FOUT2.write("chr\tstart\tstop\tfamily\t")
FOUT2.write("\t".join(["b_" + x for x in nicenames]))
FOUT2.write("\t")
FOUT2.write("\t".join(["s_" + x for x in nicenames]))
FOUT2.write("\n")

chromosomes = []
this_chromosome = None
group_count = 0
while True:
    # find peak to work on

    # check if we're still on the current chromosome
    next_chromosomes = [p.peek.seqid for p in parsers if p.peek]
    if len(next_chromosomes) == 0:
        # end of file:
        break

    if not this_chromosome in next_chromosomes:
        assert(len(set(next_chromosomes)) == 1)
        this_chromosome = next_chromosomes[0]

    # find the lowest coordinate - and pick that as a reference
    reference_peek = None
    for i, parser in enumerate(parsers):
        bmp = parser.peek
        if not bmp:
            continue
        if bmp.seqid != this_chromosome:
            # not on this chromosome - we'll pick this one up later
            continue
        if not reference_peek:
            reference_peek = bmp
        elif bmp.start < reference_peek.start:
            reference_peek = bmp

    # now we have a reference peek - start collecting peeks that
    # overlap with the reference peek - and are of the same family

    def gff_type_to_fam(s):
        return s.replace('REFS.', '')\
            .replace('UNIQUE.', '')\
            .replace('NOTUNIQ.', '')\
            .replace('(', '')\
            .replace(')', '')\
            .replace('#', '')

    # current bump stats
    bump_start = reference_peek.start
    bump_end = reference_peek.end
    bump_seqid = reference_peek.seqid
    bump_type = gff_type_to_fam(reference_peek.type)

    bump_group = []

    while True:
        found_overlap = False
        for i, parser in enumerate(parsers):
            check_bmp = parser.peek

            if not check_bmp:
                continue

            check_type = gff_type_to_fam(check_bmp.type)

            if check_bmp.seqid != bump_seqid:
                # different chromosome
                continue
            if check_bmp.start > bump_end:
                # past the current bump - ignore
                continue
            if check_type != bump_type:
                # different type (family) - ignore
                continue
            # found an overlap - same chrom, same type:
            new_bump = parser.next()
            bump_group.append(new_bump)
            bump_end = max(bump_end, new_bump.end)
            found_overlap = True
        if not found_overlap:
            break

    if not bump_group:
        # assume end of file
        break

    tags = sorted(list(set([b.tag for b in bump_group])))
    scores = collections.defaultdict(int)
    for b in bump_group:
        scores[b.tag] = max(scores[b.tag], b.score)

    # print tags
    # print scores

    ld = [{True: '1', False: '0'}[nn in tags] for nn in nicenames]
    ls = [str(scores[nn]) for nn in nicenames]

    # print ld # print ls
    if args.only_differential:
        if len(set(ld)) < 2:
            continue

    group_count += 1
    FOUT1.write("\t".join(map(str, [bump_seqid, bump_start, bump_end, bump_type, str(len(tags))] + tags)))
    FOUT1.write("\n")
    FOUT2.write("bg%08d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (group_count, bump_seqid, bump_start, bump_end, bump_type, "\t".join(ld), "\t".join(ls)))
    gff_attrs = ['ID=REFS.REGION.%s.%s.%s' % (
                    bump_type, bump_seqid, bump_start)]
    for t in tags:
        gff_attrs.append('Present_in=%s' % t)

    COREGFF.write("%s\n" % "\t".join(map(str, [
                bump_seqid, 'REFS', 'REFS.REGION.%s' % bump_type, bump_start, bump_end,
                '.', '.', '.', ";".join(gff_attrs)])))
