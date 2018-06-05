import unittest, hts, genoiser

var b: Bam
open(b, "tests/x.bam", index=true)
var depths = Fun(f:depthfun)
var funs = @[depths]

suite "genoiser-suite":
  test "test that region is set":
    check b != nil
    check b.genoiser(funs, "1", 1000000, 2000000)
    var found = false
    for v in depths.values:
      if v > 0:
        found = true
        break
    check found
    check depths.values.len == 2000000 - 1000000

  test "test that empty returns false":
    check: not b.genoiser(funs, "10", 1000000, 2000000)
