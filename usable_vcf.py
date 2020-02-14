#!/usr/bin/env python
"""
Runs commands and VCF parsing libraries to check if a vcf is usable or not
"""
import os
import sys
import logging
import argparse

import pysam
import vcf

from acebinf import cmd_exe, setup_logging

def pcmd_exe(cmd):
    """
    Wraps a cmd_exe with set -o pipefail
    """
    return cmd_exe("set -o pipefail; " + cmd)

def parse_args(args):
    parser = argparse.ArgumentParser(prog="usable_vcf", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf",
                        help="vcf to parse")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Output errors returned by commands/libraries")
    parser.add_argument("-m", "--max-entries", default=10000,
                        help="Maximum number of entries to parse (%(default)s)")
    return parser.parse_args(args)

 
def main(args):
    """
    """
    args = parse_args(args)
    setup_logging(args.verbose)
    vfn = args.vcf
    if not os.path.exists(vfn):
        logging.error("%s does not exist.")
        
    v = pysam.VariantFile(vfn)
    
    all_infos = [x for x in v.header.info]
    all_formats = [x for x in v.header.formats]
    fstr = " ".join(["%INFO/" + x for x in all_infos]) + " [" + " ".join([ "%" + x for x in all_formats]) + "]"

    failures = 0
    tests = 0
    ret = pcmd_exe("bcftools query -f '%s' %s | head -n%d" % (fstr, vfn, args.max_entries))
    tests += 1
    if ret.ret_code != 0 or ret.stderr != b"":
        logging.error("===== bcftools query failed =====")
        logging.debug(ret.stderr.decode())
        failures += 1


    # Need a way to generate/check regions (maybe multiple of them)...
    ret = pcmd_exe("bcftools view %s | head -n%d" % (vfn, args.max_entries))
    tests += 1
    if ret.ret_code != 0 or ret.stderr != b"":
        logging.error("===== bcftools view failed =====")
        logging.debug(ret.stderr.decode())
        failures += 1

    ret = pcmd_exe("vcf-validator -u -d %s" % (vfn))
    tests += 1
    if ret.ret_code != 0 or ret.stdout != "":
        logging.error("===== vcf-validator failed =====")
        logging.debug(ret.stdout)
        failures += 1


    # pysam read checker
    i = 0
    tests += 1
    try:
        for entry in v:
            i += 1
            if i > args.max_entries:
                break
            for key in entry.info:
                if key not in all_infos:
                    raise Exception("problem %s" % key)
            for sample in entry.samples:
                for key in entry.samples[sample]:
                    if key not in all_formats:
                        raise Exception("problem %s" % key)
    except Exception as e:
        logging.error("===== pysam failed ===")
        logging.debug(e)
        failures += 1

    v = vcf.VCFReader(filename=vfn)
    i = 0
    tests += 1
    try:
        for entry in v:
            i += 1
            if i > args.max_entries:
                break
            for key in entry.INFO:
                if key not in all_infos:
                    raise Exception("problem %s" % key)
                    break
            for sample in entry.samples:
                for key in sample.data._fields:
                    if key not in all_formats:
                        raise Exception("problem %s" % key)
                        break
    except Exception as e:
        logging.error("===== pyvcf failed ===")
        logging.debug(e)
        failures += 1


    logging.info("===== %d/%d passed (%d failures) for %s =====" % (tests - failures, tests, failures, vfn))
        
if __name__ == '__main__':
    main(sys.argv[1:])
