# Package

version       = "0.2.5"
author        = "Brent Pedersen"
description   = "write functions, get summaries of genomic data"
license       = "MIT"


# Dependencies

requires "hts >= 0.2.3", "docopt >= 0.6.8", "kexpr"
srcDir = "src"

bin = @["genoiser"]
installExt = @["nim"]


skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r --threads:on tests/all"

#before test:
#  exec "c2nim src/hts/private/hts_concat.h"

task docs, "Builds documentation":
  mkDir("docs"/"genoiser")
  #exec "nim doc2 --verbosity:0 --hints:off -o:docs/index.html  src/hts.nim"
  for file in listfiles("src/genoiser"):
    if file.endswith("value.nim"): continue
    if splitfile(file).ext == ".nim":
      exec "nim doc2 --verbosity:0 --hints:off -o:" & "docs" /../ file.changefileext("html").split("/", 1)[1] & " " & file

