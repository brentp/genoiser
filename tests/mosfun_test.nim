import unittest, hts, mosfun

var b: Bam
open(b, "tests/x.bam", index=true)
var depths = Fun(f:depthfun)
var funs = @[depths]

suite "mosfun-suite":
  test "test that region is set":
    check b != nil
    check b.mosfun(funs, "1", 1000000, 2000000)
    var found = false
    for v in depths.values:
      if v > 0:
        found = true
        break
    check found
    check depths.values.len == 2000000 - 1000000

  test "test that empty returns false":
    check: not b.mosfun(funs, "10", 1000000, 2000000)
