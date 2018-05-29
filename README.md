#### mosfun: write (nim) functions to get the most from your BAMs/CRAMs

[mosdepth](https://github.com/brentp/mosdepth) uses chromosome-sized arrays of
int32's to track sequencing depth. This is [fast and flexible](https://brentp.github.io/post/arrays/).

given `mosdepth` as a special-case for *depth*, `mosfun` is a general case for user-defined functions.
`mosdepth` could be implemented with `mosfun`.

An added benefit is the reduction of memory; `mosdepth` allocates an int32 array the size of each
chromosome--meaning about 1GB of memory for chromosome 1. `mosfun` can use smaller-sized chunks to
tile across each chromosome. This is important because it uses 1 array for each user-defined function.
It defaults to 8 megabase chunks as that is the smallest size with no noticeable effect on performance
in our tests. Chunk sizes down to 100KB have some, but minor effect on performance.

The idea is of `mosfun` is that it will handle all accounting, a user
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


## Speed

for maximum speed, compile with `nim c -d:release --passC:-flto --passL:-s --gc:markAndSweep src/mosfun.nim`

## CLI

The command-line interface allows running pre-specified noise filters in 2 steps. The first steps calculates the "noise" in each sample.

```
mosfun per-sample mosfun per-sample --fasta $reference results/$sample /path/to/$sample.bam # or cram
```

This will use a single thread and it will take about 1 hour for 30X bams and a bit more for crams. It will output 5 files per sample.
at the specificied prefix, in this case, `results/$sample`

After all samples have been run, then the user can `aggregate` the signal across samples. The recommended commands are:

```
chroms="$(seq 1 22) X Y"

# count number of samples at each site where 10 or more reads read with 4 or more mismatches overlapped.
echo $chroms | tr ' ' '\n' \
    | gargs -d -v -p 22 "mosfun aggregate -e 'value > 10' results/*.mosfun.{}.mismatches.bed > mismatches.{}.bed"

# count number of samples at each site where more than 1 and less than 15% of reads were interchromosomal.
echo $chroms | tr ' ' '\n' \
    | gargs -d -v -p 22 "mosfun aggregate -e '(value / depth) < 0.15 & (value > 1)' results/*.mosfun.{}.interchromosomal.bed > interchromosomal.{}.bed"

# count number of samples at each site where more than 1 and less than 15% of reads were weird.
echo $chroms | tr ' ' '\n' \
    | gargs -d -v -p 22 "mosfun aggregate -e '(value / depth) < 0.15 & (value > 1)' results/*.mosfun.{}.weird.bed > weird.{}.bed"

# count number of samples at each site where more than 1 and less than 15% of reads were soft-clipped.
echo $chroms | tr ' ' '\n' \
    | gargs -d -v -p 22 "mosfun aggregate -e '(value / depth) < 0.15 & (value > 1)' results/*.mosfun.{}.soft.bed > soft.{}.bed"

# count number of samples at each site where more than 2 where soft-clipped.
echo $chroms | tr ' ' '\n' \
    | gargs -d -v -p 22 "mosfun aggregate -e '(value > 2)' results/*.mosfun.{}.soft.bed > high-soft.{}.bed"

```
where `gargs` is available as a static binary from [here](https://github.com/brentp/gargs/releases)

This will output 1 file per chromosome, per metric. The resulting files are the desired output.


