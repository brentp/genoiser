import hts
import math
import os
import tables
import docopt
import sequtils
import strutils
import algorithm
import binaryheap
import kexpr

type 
  mrange* = tuple[start:int, stop:int, count: int]
  mchrom* = ref object
      chrom:string
      start:int
      stop:int
      count: int
  mpair* = tuple[start:int, stop:int, depth: int, value:int]
  crange* = ref object
    chrom:string
    start:int
    stop:int
    ## index of the file
    idx:int
  cend* = ref object
    chrom:string
    pos:int
    ## up (1) or down(-1)
    upd: int
    ## index of the file

proc accumulater[T](c: var seq[T]) =
  # convert from an array of start/end inc/decs to actual coverage.
  var tracker = T(0)
  for i, v in pairs(c):
    tracker += v
    c[i] = tracker

iterator ranges*[T](counts: var seq[T], chrom:string): mchrom {.inline.} =
  var last_count = counts[0]
  var last_i = 0
  for i, c in counts:
    if last_count != c:
      if last_count != 0:
        yield mchrom(chrom:chrom, start:last_i, stop:i, count:last_count.int)
      last_i = i
      last_count = c
  if last_count != 0:
    yield mchrom(chrom:chrom, start:last_i, stop:counts.len, count:last_count)

iterator mranges*[T](depth: var seq[T], values: var seq[T]): mpair {.inline.} =
  ## merge intervals+values where consecutive values and depths are unchanged
  if len(depth) != len(values):
    raise newException(ValueError, "mosfun:expected equal length arrays")

  var last_pair = (depth[0], values[0])
  var last_i = 0
  for i, d in depth:
    var v = values[i]
    if last_pair[0] != d or last_pair[1] != v:
      if last_pair[1] != 0:
        yield (last_i, i, last_pair[0].int, last_pair[1].int)
      last_i = i
      last_pair = (d, v)
  if last_pair[1] != 0:
    yield (last_i, len(depth), last_pair[0].int, last_pair[1].int)

type
  Fun* = ref object
    values*: seq[int32]
    f*: proc(aln:Record, posns:var seq[mrange])

proc mosfun*(bam: Bam, funs: seq[Fun], chrom: string, start:int, stop:int): bool =
  ## for the chromosome given, call each f in fs if skip_fun returns false and return
  ## an array for each function in fs.
  result = false

  var tid: int = -1
  var tlen = -1
  for t in bam.hdr.targets:
    if t.name == chrom:
      tid = t.tid
      tlen = int(t.length)
  if tid == -1:
    raise newException(KeyError, "chromosome not found:" & chrom)

  for i, f in funs:
    if f.values == nil or f.values.len != stop - start + 1:
      echo "creating new seq"
      f.values = new_seq[int32](stop - start)

  var posns = new_seq_of_cap[mrange](200)

  for record in bam.queryi(tid, max(0, start - 1), stop + 1):
    for i, f in funs:
      if posns.len != 0: posns.set_len(0)
      f.f(record, posns)
      for se in posns:
        var se_start = max(0, se.start - start)
        if se_start >= f.values.len:
          continue
        var se_stop = se.stop - start
        if se_stop < 0:
          continue
        if se_stop >= f.values.len:
          se_stop = f.values.len - 1
        result = true
        f.values[se_start] += int32(se.count)
        f.values[se_stop] -= int32(se.count)

  if result:
    for f in funs:
      accumulater(f.values)

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

proc eventfun*(aln:Record, posns:var seq[mrange]) =
  ## eventfun is an example of a `fun` that can be sent to `mosfun`.
  ## it sets positions where there are soft-clips, hard-clips, insertions, or deletions.
  var f = aln.flag
  if f.unmapped or f.secondary or f.supplementary or f.qcfail or f.dup: return
  var cig = aln.cigar
  if cig.len == 1: return
  var pos = aln.start

  for op in cig:
    case op.op:
    of CigarOp.soft_clip, CigarOp.hard_clip, CigarOp.insert:
      # for this function, we want the exact break-points, not the span of the event,
      # so we increment the position and the one that follows it.
      posns.add((pos, pos+1, 1))
    of CigarOp.deletion:
      posns.add((pos, pos+1, 1))
      posns.add((pos+op.len, pos+op.len+1, 1))
    else:
      discard
    if op.consumes.reference:
      pos += op.len

proc read_line(hf: ptr htsFile, kstr: ptr kstring_t, check_ok: proc(depth:int, value:int):bool, idx:int): crange {.inline.} =
  
  var iv = crange(idx:idx)
  while true:
    if hts_getline(hf, cint(10), kstr) <= 0:
      return nil

    var
      i = 0
      depth:int
      value:int
    # using an iterator of the split and setting each attribute is faster
    # than creating a slice.
    for tok in ($kstr.s).strip().split("\t"):
      if i == 0:
        iv.chrom = tok
      elif i == 1:
        iv.start = parseInt(tok)
      elif i == 2:
        iv.stop = parseInt(tok)
      elif i == 3:
        depth = parseInt(tok)
      elif i == 4:
        value = parseInt(tok)
      i+=1

    if check_ok(depth, value):
      return iv
  return nil

type sample_count {.shallow.} = ref object
  arr:seq[int8]
  chrom:string


proc sample_binary(f: string, sample_checker: proc(depth:int, value:int): bool, fai:Fai): sample_count =

  var kstr = kstring_t(l:0, m:0, s:nil)

  var fh = hts_open(f.cstring, "r")
  var v = read_line(fh, kstr.addr, sample_checker, 0)
  result = sample_count(chrom:v.chrom, arr:new_seq[int8](fai.chrom_len(v.chrom) + 1))
  while v != nil:
    when defined(debug):
      if v.stop >= result.arr.len:
        quit "got stop > chromosome length"
    result.arr[v.start].inc
    result.arr[v.stop].dec
    v = read_line(fh, kstr.addr, sample_checker, 0)
  discard hts_close(fh)


iterator aggregator*(fns: seq[string], sample_checker: func(depth:int, value:int): bool, fai:Fai): mchrom =
  var a: seq[int32]
  var chrom: string = ""

  for i, f in fns:
    if i > 0 and i mod 100 == 0:
      stderr.write_line "on file " & $i & " of " & $fns.len
    var t = sample_binary(f, sample_checker, fai)
    shallow(t.arr)
    if chrom != "" and chrom != t.chrom:
      quit "got different chroms"
    if chrom == "":
      chrom = t.chrom
      a = newSeq[int32](fai.chrom_len(chrom) + 1)
    else:
      if t.arr.len != a.len:
        quit "different length chromosomes. only send in files from same chrom."
    for i, v in t.arr:
      # we know the values should not overlap.
      when defined(debug):
        if not (v == -1 or v == 0 or v == 1): quit "got unexpected value at:" & $i & " of " & $v
      if v != 0:
        a[i] += v.int32

  accumulater(a)
  for m in ranges(a, chrom):
    yield m

proc refposns(aln:Record, posns:var seq[mrange]) {.inline.} =
  var c = aln.cigar
  # generate start, end pairs given a cigar string and a position offset.
  var pos = aln.start
  for op in c:
    var c = op.consumes
    if not c.reference:
      continue
    var olen = op.len
    if c.query:
      if len(posns) == 0 or pos != posns[len(posns)-1].stop:
        # for this function, we want the span of the event so we increment start .. stop.
        posns.add((pos, pos+olen, 1))
      else:
        posns[len(posns)-1].stop = pos + olen
    pos += olen

proc depthfun*(aln:Record, posns:var seq[mrange]) =
  ## depthfun is an example of a `fun` that can be sent to `mosfun`.
  ## it sets reports the depth at each position
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  refposns(aln, posns)

proc concordant(aln:Record): bool {.inline.} =
  if aln.tid != aln.mate_tid: return false
  # TODO: make these data-driven
  if aln.isize.abs > 600: return
  if aln.isize.abs < 50: return
  var f = aln.flag
  # check we have +- orientation.
  return f.reverse != f.mate_reverse and aln.start > aln.mate_pos == f.reverse

proc mismatchfun(aln:Record, posns: var seq[mrange]) =
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  var nm = tag[int](aln, "NM")
  if nm.isNone or nm.get < 3: return
  posns.add((aln.start, aln.stop, 1))


proc fragfun*(aln:Record, posns:var seq[mrange]) =
  ## if true proper pair, then increment the entire fragment. otherwise, increment
  ## the start and end of each read separately.
  #if aln.mapping_quality < 5: return
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  #if aln.tid != aln.mate_tid: return
  #if not f.proper_pair: return
  if aln.concordant:
    # only count the fragment once.
    if f.reverse: return
    if aln.isize < 0:
      quit "BAD"
    posns.add((aln.start, aln.start + aln.isize, 1))
  #else:
    #echo aln.flag
    #echo "disconcordant"
    #posns.add((aln.start, aln.stop, 1))


proc weird*(aln:Record, posns:var seq[mrange]) =
  ## weird increments from read-start to end for paired reads that are not in the usual orientation.
  var f = aln.flag
  if f.mate_unmapped or f.unmapped or f.secondary or f.qcfail or f.dup: return
  #if aln.stop > aln.mate_pos: return
  if aln.tid != aln.mate_tid: return
  if f.reverse == f.mate_reverse or aln.start < aln.mate_pos == f.reverse:
    posns.add((aln.start, aln.stop, 1))

proc mq0fun*(aln:Record, posns:var seq[mrange]) =
  ## this is an example function that increments all reference locations with mapping-quality 0.
  if aln.mapping_quality != 0: return
  var f = aln.flag
  if f.unmapped or f.qcfail or f.dup: return
  refposns(aln, posns)

proc interchromosomal_splitter*(aln:Record, posns:var seq[mrange]) =
  for sp in aln.splitters:
    if sp.qual == 0: continue
    if sp.chrom == aln.chrom: continue
    softfun(aln, posns)
    return

proc interchromosomal*(aln:Record, posns:var seq[mrange]) =
  if aln.mapping_quality == 0: return
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  if aln.b.core.tid == aln.b.core.mtid:
    if abs(aln.start - aln.mate_pos) > 10000000:
      refposns(aln, posns)
    return
  if aln.b.core.tid == -1 or aln.b.core.mtid == -1: return
  # on different chroms
  refposns(aln, posns)

proc aggregate_main(argv: seq[string]) =

  let doc = format("""

  Usage: aggregate [options] <per-sample-output>...

Arguments:

  <per-sample-output>          the txt output from mosfun per-sample

Options:
  -f --fasta <reference>        indexed fasta reference file required to get chromosome lengths.
  -e --expr <string>            per-sample filtering expression [default: "(value > 0) & ((value / depth) > 0.1)"]
  -h --help                     show help
  """)

  let args = docopt(doc, argv=argv)

  var expr:string = ($args["--expr"]).strip(chars={'"'})

  var fai:Fai
  if $args["--fasta"] == "nil":
    quit "--fasta is required"
  if not open(fai, $args["--fasta"]):
    quit "invalid --fasta: " & $args["--fasta"]

  var ex = expression(expr)
  if ex.error() != 0:
    quit "error with expression"
  var vc = "value".cstring
  var dc = "depth".cstring

  proc sample_ok(depth:int, value:int): bool =
    discard ke_set_int(ex.ke, vc, value)
    discard ke_set_int(ex.ke, dc, depth)
    result = ex.get_bool()
    if ex.error() != 0:
      quit format("expresion error with value:$# depth:$#", value, depth)


  for iv in aggregator(@(args["<per-sample-output>"]), sample_ok, fai):
    if iv.count == 0: continue
    stdout.write_line iv.chrom & "\t" & intToStr(iv.start) & "\t" & intToStr(iv.stop) & "\t" & intToStr(iv.count)

proc writefn(a:Fun, depths:Fun, fh:File, chrom:string, start:int, min_depth:int=0, min_value:int=1) =
  for m in mranges(depths.values, a.values):
    if m.depth < min_depth: continue
    if m.value < min_value: continue
    fh.write_line chrom & "\t" & intToStr(m.start + start) & "\t" & intToStr(m.stop + start) & "\t" & intToStr(m.depth) & "\t" & intToStr(m.value)

proc myopen(path:string): File =
  var fh:File
  if not open(fh, path, fmWrite):
    quit(2)
  return fh

proc per_sample_main(argv: seq[string]) =

  let doc = format("""

  Usage: per-sample [options] <PREFIX> <BAM-or-CRAM>

Arguments:

  <PREFIX>       output path prefix
  <BAM-or-CRAM>  the alignment file.

Options:

  -f --fasta     <path>           fasta file for use with CRAM files.
  -d --min-depth <int>            min depth to report a region [default: 5]
  -v --min-value <int>            min count to report a region [default: 1]
  -h --help                       show help
  """)

  let args = docopt(doc, argv=argv)

  var fasta: cstring = nil
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])


  let prefix = $args["<PREFIX>"]
  let fbam = $args["<BAM-or-CRAM>"]
  let min_depth = parseInt($args["--min-depth"])
  let min_value = parseInt($args["--min-value"])
  var b:Bam
  open(b, fbam, index=true, fai=fasta)
  if b == nil:
    quit "coudn't open bam: " & fbam
  if b.idx == nil:
    quit "coudn't open bam index for " & fbam

  var L = 5000000
  var depths = Fun(values:new_seq[int32](L+1), f:depthfun)
  var softs = Fun(values:new_seq[int32](L+1), f:softfun)
  var mq0 = Fun(values:new_seq[int32](L+1), f:mq0fun)
  var weirds = Fun(values:new_seq[int32](L+1), f:weird)
  var misms = Fun(values:new_seq[int32](L+1), f:mismatchfun)
  var inters = Fun(values:new_seq[int32](L+1), f:interchromosomal)
  var fns = @[depths, softs, mq0, weirds, misms, inters]

  for target in b.hdr.targets:

    var fhs: seq[File]

    stderr.write_line target.name
    var start = 0
    while start < target.length.int:
      var stop = min(start + L, target.length.int)
      #echo "start..stop:", $start, "..", $stop
      if stop - start + 1 != fns[0].values.len:
        if start != 0 and stop != target.length.int:
          echo "resizing"
        for f in fns:
          f.values.set_len(stop - start + 1)
          zeroMem(f.values[0].addr.pointer, f.values.len * sizeof(f.values[0]))

      if mosfun(b, fns, target.name, start, stop):
        if fhs == nil:
          fhs = @[
            # NOTE: fragile. make sure these are same orders as fns array above.
            myopen(prefix & ".mosfun." & target.name & ".soft.bed"),
            myopen(prefix & ".mosfun." & target.name & ".mq0.bed"),
            myopen(prefix & ".mosfun." & target.name & ".weird.bed"),
            myopen(prefix & ".mosfun." & target.name & ".mismatches.bed"),
            myopen(prefix & ".mosfun." & target.name & ".interchromosomal.bed"),
          ]

        writefn(softs, depths, fhs[0], target.name, start, min_depth=min_depth, min_value=min_value)
        writefn(mq0, depths, fhs[1], target.name, start, min_depth=min_depth, min_value=min_value)
        writefn(weirds, depths, fhs[2], target.name, start, min_depth=min_depth, min_value=min_value)
        writefn(misms, depths, fhs[3], target.name, start, min_depth=min_depth, min_value=min_value)
        writefn(inters, depths, fhs[4], target.name, start, min_depth=min_depth, min_value=min_value)
        for f in fns:
          zeroMem(f.values[0].addr.pointer, f.values.len * sizeof(f.values[0]))

      start += L

    for f in fhs:
      f.close()


var progs = {
   "per-sample": per_sample_main,
   "aggregate": aggregate_main
}.toTable

var helps = {
   "per-sample": "calculate various metrics for a per-sample BAM/CRAM",
   "aggregate": "aggregate many sample outputs of per-sample into single file"
}.toTable

const version = "0.0.2"

proc main() =

  var args = commandLineParams()
  if len(args) < 1 or not progs.contains(args[0]):
    var hkeys = toSeq(keys(helps))
    sort(hkeys, proc(a, b: string): int =
      if a < b: return -1
      else: return 1
      )
    echo format("\nmosfun utility programs.\nversion: $#\n", version)

    for k in hkeys:
      echo format("	• $1: $2", k & repeat(" ", 20 - len(k)), helps[k])
    echo ""
  else:
    var p = args[0]; args.delete(0)
    progs[p](args)

when isMainModule:
  main()
