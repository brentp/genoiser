# Package

version       = "0.1.0"
author        = "Brent Pedersen"
description   = "write functions, get summaries of genomic data"
license       = "MIT"


# Dependencies

requires "hts >= 0.1.9", "docopt", "binaryheap", "kexpr"
srcDir = "src"

bin = @["mosfun"]


skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"

#before test:
#  exec "c2nim src/hts/private/hts_concat.h"

task docs, "Builds documentation":
  mkDir("docs"/"mosfun")
  #exec "nim doc2 --verbosity:0 --hints:off -o:docs/index.html  src/hts.nim"
  for file in listfiles("src/mosfun"):
    if file.endswith("value.nim"): continue
    if splitfile(file).ext == ".nim":
      exec "nim doc2 --verbosity:0 --hints:off -o:" & "docs" /../ file.changefileext("html").split("/", 1)[1] & " " & file

