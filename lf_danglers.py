#!/usr/bin/env python

import sys
import subprocess as sp
import argparse
from builtins import zip


def run_bowtie2(database, fastq, preset, threads):
    command = "bowtie2 --%s --reorder --no-head --quiet -p %d -x %s -U %s" % (
        preset, threads, database, fastq)
    process = sp.Popen(command, shell=True, stdout=sp.PIPE)
    for line in process.stdout:
        if line[0] == '@':
            continue
        yield line.split()


def align_paired_reads(args):
    return zip(run_bowtie2(args.bowtie2_database,
                           args.forward_reads,
                           args.bowtie_preset,
                           args.threads),
               run_bowtie2(args.bowtie2_database,
                           args.reverse_reads,
                           args.bowtie_preset,
                           args.threads))


def identify_danglers(aligned_pairs):
    for forward, reverse in aligned_pairs:
        forward_flag = int(forward[1])
        reverse_flag = int(reverse[1])
        if forward_flag & 0x4 == reverse_flag & 0x4:
            # either both map, or both do not map
            # not dangling
            continue

        if forward_flag & 0x4:
            # forward is dangling
            tag = 'R' + {0: '+', 0x10: '-'}[reverse_flag & 0x10]
            header = ('@%s__%s__%s\n' % (reverse[2], tag, forward[0]))
            sequence = ('%s\n+\n%s\n' % (forward[9], forward[10]))
            dangler = header + sequence
            yield dangler

        else:
            # reverse is dangling
            tag = 'F' + {0: '+', 0x10: '-'}[forward_flag & 0x10]
            header = ('@%s__%s__%s\n' % (forward[2], tag, reverse[0]))
            sequence = ('%s\n+\n%s\n' % (reverse[9], reverse[10]))
            dangler = header + sequence
            yield dangler


def write_danglers(danglers, args):
    with open(args.output, 'w') as output_file:
        output_file.writelines(danglers)
        output_file.close()


def parse_args(args):
    parser = argparse.ArgumentParser(description='find danglers')
    parser.add_argument('-t', '--threads', default=10, type=int)
    parser.add_argument('bowtie2_database')
    parser.add_argument('forward_reads')
    parser.add_argument('reverse_reads')
    parser.add_argument('output')
    parser.add_argument('-p', '--bowtie_preset', default='fast')
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    aligned_pairs = align_paired_reads(args)
    danglers = identify_danglers(aligned_pairs)
    write_danglers(danglers, args)


if __name__ == '__main__':
    main()
