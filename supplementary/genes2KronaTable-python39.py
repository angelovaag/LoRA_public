#!/usr/bin/env python3
#### this script is from: https://github.com/EnvGen/metagenomics-workshop/tree/master/in-house
#### which is from: https://metagenomics-workshop.readthedocs.io/en/latest/annotation/summarize.html
## then it was converted into python3 version by ChatGPT (tested, works)
#################################

import sys
import csv
import numpy as np

def ReadLengths(f):
    d = {}
    with open(f, 'r') as hin:
        for line in hin:
            line = line.rstrip()
            g, l = line.rsplit()
            l = float(l)
            d[g] = l
    return d

def ReadLimits(f):
    limit = {}
    with open(f, 'r') as hin:
        for line in hin:
            line = line.rstrip()
            limit[line.rsplit()[-1]] = ""
    return limit

def ReadMap(f):
    d = {}
    with open(f, 'r') as hin:
        hincsv = csv.reader(hin, delimiter='\t')
        for row in hincsv:
            ann = row[1]
            parent = row[0]
            try:
                d[ann].append(parent)
            except KeyError:
                d[ann] = [parent]
    return d

def ReadCoverage(f, lengths):
    d = {}
    try:
        with open(f, 'r') as hin:
            hincsv = csv.reader(hin, delimiter='\t')
            for row in hincsv:
                g = row[0]
                try:
                    c = float(row[1])
                except ValueError:
                    continue
                if len(lengths.keys()) > 0:
                    try:
                        l = lengths[g]
                    except KeyError:
                        sys.exit("No length found for gene " + g + "\n")
                else:
                    l = 1
                c = float(c) / l
                d[g] = c
    except TypeError:
        return d
    return d

def ReadHierarchy(f):
    d = {}
    with open(f, 'r') as hin:
        hincsv = csv.reader(hin, delimiter='\t')
        l = []
        for row in hincsv:
            id = row[0]
            name = row[1]
            try:
                hiers = row[2]
            except IndexError:
                hiers = "Unknown"
            d[id] = hiers.split("|") + [name]
            l.append(len(d[id]))
    return max(l), d

def Calculate(hier_c, operation):
    hier_sdev = {}
    if operation == "sum":
        function = np.sum
    elif operation == "mean" or operation == "meanhalf":
        function = np.mean
    for hier, l in hier_c.items():
        hier_sdev[hier] = 0.0
        if operation == "meanhalf":
            l.sort()
            l = l[len(l) // 2:]
        hier_c[hier] = function(l)
        hier_sdev[hier] = np.std(l)
    return hier_sdev, hier_c

def CalcHierarchy(mapping, annotations, coverage, operation, limit, verbose):
    ann_c = {}
    if len(coverage.keys()) > 0:
        cov = True
    else:
        cov = False
    for annotation, l in annotations.items():
        covlist = []
        if cov:
            for gene in l:
                try:
                    gene_cov = coverage[gene]
                except KeyError:
                    sys.exit("ERROR: Could not find coverage for gene " + str(gene) +
                             ". Are you sure you have coverage information?\n")
                covlist.append(gene_cov)
            ann_c[annotation] = np.mean(covlist)
        else:
            ann_c[annotation] = len(l)

    hier_c = {}
    for annotation, count in ann_c.items():

        try:
            parents = mapping[annotation]
        except KeyError:
            if verbose:
                sys.stderr.write("WARNING: Could not find hierarchy parent for " + annotation + "\n")
            continue
        for parent in parents:
            if limit and not parent in limit:
                if verbose:
                    sys.stderr.write("Skipping parent " + parent + "\n")
                continue
            try:
                hier_c[parent].append(count)
            except KeyError:
                hier_c[parent] = [count]
    hier_sdev, hier_c = Calculate(hier_c, operation)
    return hier_sdev, hier_c

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str,
                        help="Tab-delimited file with gene_ids in the first column and gene_annotations in the second column")
    parser.add_argument("-m", "--mapfile", type=str,
                        help="Tab-delimited file mapping each gene_annotation to the nearest parent hierarchy level (e.g., pathway to enzyme)")
    parser.add_argument("-H", "--hierarchy", type=str,
                        help="Hierarchy file for parent levels")
    parser.add_argument("-l", "--limit", type=str,
                        help="Limit calculations to only this list of parent hierarchies. E.g., a list of predicted pathways only")
    parser.add_argument("-O", "--operation", type=str, default="mean",
                        help="Specify how to do calculations, either 'mean' (default), 'sum' or 'meanhalf' where averages will be calculated from the most abundant half of annotations")
    parser.add_argument("-o", "--outfile", type=str,
                        help="Write results to outfile. Defaults to stdout")
    parser.add_argument("-n", "--name", type=str,
                        help="OPTIONAL: Name to assign to outfile. Defaults to the name of infile")
    parser.add_argument("-c", "--coverage", type=str,
                        help="OPTIONAL: Supply a file of coverage for each gene in the infile")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Run in verbose mode")
    parser.add_argument("-s", "--singlecol", action="store_true",
                        help="Write only one annotation column for the first parent hierarchy (e.g., pathway)")
    parser.add_argument("-d", "--sdev", type=str,
                        help="Write the standard deviation of each first parent hierarchy level to this file")
    parser.add_argument("-L", "--lengthnorm", type=str,
                        help="Provide a file with lengths for genes to normalize coverage by")
    parser.add_argument("-g", "--missing", type=str, default="missing_pathways.txt",
                        help="File name for missing pathways in the hierarchy file")
    args = parser.parse_args()

    err_f = open(args.missing, "w")

    if not args.infile or not args.mapfile or not args.hierarchy:
        sys.exit(parser.print_help())

    if args.limit:
        limit = ReadLimits(args.limit)
    else:
        limit = []

    if args.lengthnorm:
        lengths = ReadLengths(args.lengthnorm)
    else:
        lengths = {}

    mapping = ReadMap(args.mapfile)
    annotations = ReadMap(args.infile)
    coverage = ReadCoverage(args.coverage, lengths)
    max_hier, hierarchy = ReadHierarchy(args.hierarchy)
    hier_sdev, hier_counts = CalcHierarchy(mapping, annotations, coverage, args.operation, limit, args.verbose)

    if args.outfile:
        hout = open(args.outfile, 'w')
    else:
        hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter='\t')

    if args.singlecol:
        out = ["X"]
    else:
        out = ["X"] + ["Level" + str(i) for i in range(1, max_hier + 1)]
    houtcsv.writerow(out)

    if args.sdev:
        sdevout = open(args.sdev, 'w')
        sdevout.write("X.sdev\t" + args.name + "\n")

    for hier, count in hier_counts.items():
        out = [count]
        try:
            h = hierarchy[hier]
        except KeyError:
            h = ["Unknown"]
        if args.singlecol:
            out = [hier]
        else:
            out += hierarchy[hier]

        try:
            sdevout.write(hier + "\t" + str(hier_sdev[hier]) + "\n")
        except NameError:
            pass
        houtcsv.writerow(out)
    hout.close()

    try:
        sdevout.close()
    except NameError:
        pass

if __name__ == "__main__":
    main()
