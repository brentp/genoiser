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

  ## fun takes a variant and appends ranges into posns that should be incremented in an array.
  fun* = proc(aln:Record, posns:var seq[mrange])

proc accumulater[T](c: var seq[T]) =
  # convert from an array of start/end inc/decs to actual coverage.
  var tracker = T(0)
  for i, v in pairs(c):
    tracker += v
    c[i] = tracker

iterator mranges*[T](depth: var seq[T], values: var seq[T]): mpair =
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
      
proc mosfun*(bam: Bam, chrom: string, counts: var seq[seq[int32]], fs: seq[fun]) =
  ## for the chromosome given, call each f in fs if skip_fun returns false and return
  ## an array for each function in fs.

  var tid: int = -1
  var tlen = -1
  for t in bam.hdr.targets:
    if t.name == chrom:
      tid = t.tid
      tlen = int(t.length)
  if tid == -1:
    raise newException(KeyError, "chromosome not found:" & chrom)

  if counts.len != fs.len:
    counts = new_seq[seq[int32]](fs.len)
  for i, c in counts:
    if c.len != tlen:
      echo "creating new seq"
      counts[i] = new_seq[int32](tlen)

  var posns = new_seq_of_cap[mrange](200)

  for record in bam.queryi(tid, 0, tlen):
    for i, f in fs:
      if posns.len != 0: posns.set_len(0)
      f(record, posns)
      for se in posns:
        counts[i][se.start] += int32(se.count)
        counts[i][se.stop] -= int32(se.count)

  for c in counts.mitems:
    accumulater(c)

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


iterator genstream(fns: seq[string], sample_checker: proc(depth:int, value:int): bool): crange =
  ## merge all files into a single sorted stream of intervals (after filtering by sample_checker).
  var kstr : kstring_t
  kstr.l = 0
  var fhs = new_seq[ptr htsFile](len(fns))
  var q = newHeap[crange]() do (a, b: crange) -> int:
    a.start - b.start

  for i, fn in fns:
    fhs[i] = hts_open(fn.cstring, "r")
    var v = read_line(fhs[i], kstr.addr, sample_checker, i)
    if v != nil:
      q.push(v)

  while q.size != 0:
    # pop the leftest interval.
    var left = q.pop
    yield left
    var v = read_line(fhs[left.idx], kstr.addr, sample_checker, left.idx)
    if v != nil:
      q.push(v)

proc merge_pop(heap: var Heap[cend]): cend {.inline.} =
  ## if pop off the heap until we get a novel position. this
  ## makes sure we don't yield multiple identical positions from
  ## gense
  result = heap.pop
  while heap.size > 0 and heap.peek.pos == result.pos:
    var tmp = heap.pop
    result.upd += tmp.upd

iterator gense(fns: seq[string], sample_checker: proc(depth:int, value:int): bool): cend =
  ## convert start, end into position with incrment of +1 for start, -1 for end.
  var heap = newHeap[cend]() do (a, b: cend) -> int:
    a.pos - b.pos

  for iv in genstream(fns, sample_checker):
    # start is +1, end is -1, can yield a start once the end is passed.
    while heap.size > 0 and heap.peek.pos < iv.start:
      yield heap.merge_pop
    heap.push(cend(chrom:iv.chrom, pos:iv.start, upd:1))

    #while heap.size > 0 and heap.peek.pos < iv.stop:
    #  yield heap.pop
    var cc = cend(chrom:iv.chrom, pos:iv.stop, upd:(-1))
    heap.push(cc)

  while heap.size > 0:
    yield heap.merge_pop

iterator aggregator0*(fns: seq[string], sample_checker: proc(depth:int, value:int): bool): mchrom =
  # 1. creating a single merged, sorted stream of intervals from all files
  # 2. push start (+1) and end (-1) onto a queue.
  # 3. for each start added to queue. we can pop all items in the heap positioned below it.
  # 4. track cumulative sum of +1 and -1
  var 
    count = 0
    last_pos = 0

  for iv in gense(fns, sample_checker):
    yield mchrom(chrom:iv.chrom, start:last_pos, stop:iv.pos, count:count)
    last_pos = iv.pos
    count += iv.upd

iterator aggregator*(fns: seq[string], sample_checker: proc(depth:int, value:int): bool): mchrom =
  var 
    last:mchrom

  for m in aggregator0(fns, sample_checker):
    # extend last interval
    if last == nil:
      last = m
      continue
    if m.count == last.count and m.start == last.stop:
      last.stop = m.stop
    else:
      yield last
      last = mchrom(chrom:m.chrom, start:m.start, stop:m.stop, count:m.count)

  yield last

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

proc fragfun*(aln:Record, posns:var seq[mrange]) =
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  if aln.mapping_quality < 20: return
  #if not f.proper_pair: return

  if aln.stop > aln.mate_pos: return
  if aln.isize > 400: return
  if aln.isize < 200: return
  posns.add((aln.start, aln.start + aln.isize, 1))
  echo aln.start, " ", aln.stop, " ", aln.mate_pos, " ", aln.start + aln.isize, " ", aln.isize, " ", f.proper_pair

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
    #if aln.b.core.tid == -1:
    #  return
    #interchromosomal_splitter(aln, posns)
    return
  if aln.b.core.tid == -1 or aln.b.core.mtid == -1: return
  refposns(aln, posns)

proc aggregate_main(argv: seq[string]) =

  let doc = format("""

  Usage: aggregate [options] <per-sample-output>...

Arguments:

  <per-sample-output>          the txt output from mosfun per-sample

Options:
  -e --expr <string>            per-sample filtering expression [default: "(value > 0) & ((value / depth) > 0.1)"]
  -h --help                     show help
  """)

  let args = docopt(doc, argv=argv)

  var expr:string = ($args["--expr"]).strip(chars={'"'})

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
  for iv in aggregator(@(args["<per-sample-output>"]), sample_ok):
    stdout.write_line iv.chrom & "\t" & intToStr(iv.start) & "\t" & intToStr(iv.stop) & "\t" & intToStr(iv.count)


proc per_sample_main(argv: seq[string]) =

  let doc = format("""

  Usage: per-sample [options] <PREFIX> <BAM-or-CRAM>

Arguments:

  <BED>          output path prefix
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

  for target in b.hdr.targets:
    if target.name != "chr17": continue
    var softs = new_seq[int32](target.length); shallow(softs)
    var mq0 = new_seq[int32](target.length); shallow(mq0)
    var depths = new_seq[int32](target.length); shallow(depths)
    var inters = new_seq[int32](target.length); shallow(inters)
    var splits = new_seq[int32](target.length); shallow(splits)
    var frag_depths = new_seq[int32](target.length); shallow(frag_depths)

    var counts = @[softs, mq0, depths, inters, splits, frag_depths]
    mosfun(b, target.name, counts, @[fun(softfun), fun(mq0fun), fun(depthfun), fun(interchromosomal), fun(interchromosomal_splitter), fun(fragfun)])

    var fs, f0, fi, fspl, ffrag, fdepth:File
    if not open(fs, prefix & "." & target.name & ".soft.bed", fmWrite):
      quit(2)
    if not open(f0, prefix & "." & target.name & ".mq0.bed", fmWrite):
      quit(2)
    if not open(fdepth, prefix & "." & target.name & ".depth.bed", fmWrite):
      quit(2)
    if not open(fi, prefix & "." & target.name & ".interchrom.bed", fmWrite):
      quit(2)
    if not open(fspl, prefix & "." & target.name & ".splits.bed", fmWrite):
      quit(2)
    if not open(ffrag, prefix & "." & target.name & ".frag.bed", fmWrite):
      quit(2)


    for m in mranges(depths, mq0):
      if m.depth < min_depth: continue
      if m.value < min_value: continue
      f0.write_line target.name & "\t" & intToStr(m.start) & "\t" & intToStr(m.stop) & "\t" & intToStr(m.depth) & "\t" & intToStr(m.value)

    for m in mranges(depths, splits):
      if m.depth < min_depth: continue
      if m.value < min_value: continue
      fspl.write_line target.name & "\t" & intToStr(m.start) & "\t" & intToStr(m.stop) & "\t" & intToStr(m.depth) & "\t" & intToStr(m.value)

    for m in mranges(depths, softs):
      if m.depth < min_depth: continue
      if m.value < min_value: continue
      fs.write_line target.name & "\t" & intToStr(m.start) & "\t" & intToStr(m.stop) & "\t" & intToStr(m.depth) & "\t" & intToStr(m.value)

    for m in mranges(depths, inters):
      if m.depth < min_depth: continue
      if m.value < min_value: continue
      fi.write_line target.name & "\t" & intToStr(m.start) & "\t" & intToStr(m.stop) & "\t" & intToStr(m.depth) & "\t" & intToStr(m.value)

    for m in mranges(depths, frag_depths):
      #if m.depth < min_depth: continue
      #if m.value < min_value: continue
      ffrag.write_line target.name & "\t" & intToStr(m.start) & "\t" & intToStr(m.stop) & "\t" & intToStr(m.depth) & "\t" & intToStr(m.value)

    GC_fullCollect()
    fi.close(); f0.close(); fs.close(); fspl.close(); ffrag.close()
    #echo "debug: breaking early"
    #break


var progs = {
   "per-sample": per_sample_main,
   "aggregate": aggregate_main
}.toTable

var helps = {
   "per-sample": "calculate various metrics for a per-sample BAM/CRAM",
   "aggregate": "aggregate many sample outputs of per-sample into single file"
}.toTable

const version = "0.0.1"

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
