#### mosfun: write (nim) functions to get the most from your BAMs/CRAMs

[mosdepth](https://github.com/brentp/mosdepth) uses chromosome-sized arrays of
int32's to track sequencing depth. This is [fast and flexible](https://brentp.github.io/post/arrays/).

given `mosdepth` as a special-case for *depth*, `mosfun` is a general case for user-defined functions.
`mosdepth` could be implemented with `mosfun`.

This library is in progress. The idea is that `mosfun` handles all accounting, a user
simply defines a [nim](https://nim-lang.org) function that takes an alignment and then
indicates which genomic positions to increment. For example, to calculate depth, this user
function would increment from start to end:

```Nim
proc depthfun*(aln:Record, posns:var seq[mrange]) =
  ## depthfun is an example of a `fun` that can be sent to `mosfun`.
  ## it increments from aln.start to aln.stop of passing reads.
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  posns.add((aln.start, aln.stop, 1))
```

The `posns` value is sent to the function by `mosfun` and the user-defined function
can add to it as many elements as desired. In this case it increments from `aln.start`
to `aln.stop` by `1`. It can inrement by any integer value.

The user could also choose to increment any soft or hard-clip location:

```Nim
proc softfun*(aln:Record, posns:var seq[mrange]) =
  ## softfun an example of a `fun` that can be sent to `mosfun`.
  ## it sets positions where there are soft-clips
  var f = aln.flag
  if f.unmapped or f.secondary or f.supplementary or f.qcfail or f.dup: return
  var cig = aln.cigar
  if cig.len == 1: return
  var pos = aln.start

  for op in cig:
    if op.op == CigarOp.soft_clip or op.op == CigarOp.hard_clip:
      # for this function, we want the exact break-points, not the span of the event,
      # so we increment the position and the one that follows it.
      posns.add((pos, pos+1, 1))
    if op.consumes.reference:
      pos += op.len
```

### Utility

This library provides the machinery. Other command-line tools will use this for more obviously useful things.

