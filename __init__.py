
import wrappers
import statistics
import fits
import schedulers

import helpers
import IO


def create_wrapper(codestring, execline, optsdict=None):

  if codestring == "TRI3DST":
    return wrappers.TRI3DST.TRI3DST_Wrapper(execline)

  if codestring == "SDTRIMSP":
    return wrappers.SDTRIMSP.SDTRIMSP_Wrapper(execline)

  if codestring == "PARCAS":
    return wrappers.PARCAS.PARCAS_Wrapper(execline)
