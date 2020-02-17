#!/usr/bin/env python
"""
Runs commands and VCF parsing libraries to check if a vcf is usable or not
"""
import os
import sys
import time
import signal
import logging
import argparse
import datetime
import tempfile
import warnings
import subprocess

from collections import namedtuple

import vcf
import pysam


def setup_logging(debug=False, stream=sys.stderr, log_format=None):
    """
    Default logger
    """
    logLevel = logging.DEBUG if debug else logging.INFO
    if log_format is None:
        log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=stream, level=logLevel, format=log_format)
    logging.info("Running %s", " ".join(sys.argv))

    # pylint:disable=unused-argument
    def sendWarningsToLog(message, category, filename, lineno, *args, **kwargs):
        """
        Put warnings into logger
        """
        logging.warning('%s:%s: %s:%s', filename, lineno, category.__name__, message)
        return
    warnings.showwarning = sendWarningsToLog

class Alarm(Exception):

    """ Alarm Class for command timeouts """
    pass


def alarm_handler(signum, frame=None):  # pylint: disable=unused-argument
    """ Alarm handler for command timeouts """
    raise Alarm

cmd_result = namedtuple("cmd_result", "ret_code stdout stderr run_time")


def cmd_exe(cmd, timeout=-1):
    """
    Executes a command through the shell.
    timeout in minutes! so 1440 mean is 24 hours.
    -1 means never
    returns namedtuple(ret_code, stdout, stderr, datetime)
    where ret_code is the exit code for the command executed
    stdout/err is the Standard Output Error from the command
    and datetime is a datetime object of the execution time
    """
    t_start = time.time()
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, close_fds=True,
                            preexec_fn=os.setsid)
    signal.signal(signal.SIGALRM, alarm_handler)
    if timeout > 0:
        signal.alarm(int(timeout * 60))
    try:
        stdoutVal, stderrVal = proc.communicate()
        signal.alarm(0)  # reset the alarm
    except Alarm:
        os.killpg(proc.pid, signal.SIGTERM)
        proc.kill()
    	t_end = time.time()
        stdoutVal, stderrVal = proc.communicate()
    	stdoutVal = bytes.decode(stdoutVal)
    	ret = cmd_result(214, stdoutVal, stderrVal, datetime.timedelta(seconds=int(t_end - t_start)))
        return ret
    t_end = time.time()

    stdoutVal = bytes.decode(stdoutVal)
    retCode = proc.returncode
    ret = cmd_result(retCode, stdoutVal, stderrVal, datetime.timedelta(seconds=int(t_end - t_start)))
    return ret

def pcmd_exe(cmd):
    """
    Wraps a cmd_exe with set -o pipefail
    """
    out = tempfile.NamedTemporaryFile(mode='w')
    out.write("set -o pipefail; " + cmd)
    return cmd_exe("bash " + out.name)

def parse_args(args):
    parser = argparse.ArgumentParser(prog="usable_vcf", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf",
                        help="vcf to parse")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Output errors returned by commands/libraries")
    parser.add_argument("-m", "--max-entries", default=100, type=int,
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

    # header doesn't count as entries
    max_entries = args.max_entries + str(v.header).count('\n')
    
    failures = 0
    tests = 0
    ret = pcmd_exe("bcftools query -f '%s\n' %s | head -n%d" % (fstr, vfn, max_entries))
    tests += 1
    # 141 is from the 'head'
    if (ret.ret_code != 0 and ret.ret_code != 141) or ret.stderr != b"":
        logging.error("===== bcftools query failed =====")
        logging.debug(ret.stderr.decode())
        failures += 1


    # Need a way to generate/check regions (maybe multiple of them)...
    ret = pcmd_exe("bcftools view %s | head -n%d" % (vfn, max_entries))
    tests += 1
    if (ret.ret_code != 0 and ret.ret_code != 141) or ret.stderr != b"":
        logging.error("===== bcftools view failed =====")
        logging.debug(ret.stderr.decode())
        failures += 1

    # no need to pipefail here
    ret = cmd_exe("vcf-validator -d %s" % (vfn), timeout=args.max_entries/60.0)
    tests += 1
    if (ret[0] != 0 and ret[0] != 141) or ret[2] != "":
        logging.error("===== vcf-validator failed =====")
        logging.debug(ret[2])
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
            samp_cnt = 0
            for sample in entry.samples:
                if samp_cnt >= args.max_entries: continue
                samp_cnt += 1
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
            samp_cnt = 0
            for sample in entry.samples:
                if samp_cnt >= args.max_entries: continue
                samp_cnt += 1
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
