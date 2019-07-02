{.experimental.}

import hts
import math
import os
import tables
import docopt
import sequtils
import strutils
import algorithm
import threadpool
import kexpr

type
  mrange* = tuple[start:int, stop:int, count: int]
  mchrom* = ref object
      chrom:string
      start:int
      stop:int
      count: int
  mpair* = tuple[start:int, stop:int, depth: int, value:int]
  crange* = object
    chrom:string
    start:int
    stop:int
    ## index of the file
    idx:int

proc cumulative_sum*[T](c: var seq[T]) =
  # convert from an array of start/end inc/decs to actual coverage.
  # tracker allows for overflow/underflow for the data type T and just
  # keeps the vlaue at T.high until it drops back down again.
  var tracker:int64
  for i, v in pairs(c):
    tracker += v.int64
    if tracker > T.high.int64:
      c[i] = T.high
      #tot = tracker.float64
    elif tracker < T.low.int64:
      c[i] = T.low.abs
    else:
      c[i] = T(tracker.abs)

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
    yield mchrom(chrom:chrom, start:last_i, stop:counts.len, count:last_count.int)

iterator mranges*[T](depth: var seq[T], values: var seq[T]): mpair {.inline.} =
  ## merge intervals+values where consecutive values and depths are unchanged
  if len(depth) != len(values):
    raise newException(ValueError, "genoiser:expected equal length arrays")

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
  Fun* {.shallow.} [T]  = ref object
    values*: seq[T]
    f*: proc(aln:Record, posns:var seq[mrange])


proc genoiser*[T](bam_path:string, reference: string, funcs: seq[Fun[T]], chrom: string, threads:int): bool =
    var bam:Bam
    open(bam, bam_path, fai=reference, index=true)
    if bam == nil:
        quit "couldn't open bam path"
    var target:Target
    for tgt in bam.hdr.targets:
        if tgt.name == chrom:
            target = tgt
    if target == nil:
        raise newException(ValueError, "chromosome not found:" & chrom)
    for i, f in funcs:
      if f.values.len != target.length.int + 1:
        f.values = newSeq[T](target.length.int + 1)


    var channel : Channel[seq[Fun[T]]]

    proc consumer(funcs: seq[Fun[T]]) =
        echo funcs.len

    channel.open()



proc genoiser*[T](bam: Bam, funs: seq[Fun[T]], chrom: string, start:int, stop:int, cumsum:bool=true): bool =
  ## for the chromosome given, call each f in fs if skip_fun returns false and return
  ## an array for each function in fs.
  ## if cumsum is false, the array of inc/decrement values will be returned, not the cumsum. this
  ## allows combining multiple passes.
  result = false

  var tid: int = -1
  var tlen = -1
  for t in bam.hdr.targets:
    if t.name == chrom:
      tid = t.tid
      tlen = int(t.length)
  if tid == -1:
    raise newException(KeyError, "chromosome not found:" & chrom)

  var posns = new_seq_of_cap[mrange](200)

  for record in bam.query(tid, max(0, start - 1), stop + 1):
    for i, f in funs:
      if posns.len != 0: posns.set_len(0)
      f.f(record, posns)
      for se in posns:
        if f.values.len != stop - start + 1:
          echo "genoiser:creating new seq"
          f.values = new_seq[T](stop - start)

        var se_start = max(0, se.start - start)
        if se_start >= f.values.len:
          continue
        var se_stop = se.stop - start
        if se_stop < 0:
          continue
        if se_stop >= f.values.len:
          se_stop = f.values.len - 1
        result = true

        var
          sv = f.values[se_start].int64
          ev = f.values[se_stop].int64
          c = se.count

        if sv + se.count.int64 > T.high.int64:
          c = int(T.high.int64 - sv)

        if ev - c.int64 < T.low.int64:
          c = int(T.low.int64 - ev)

        if ev - c.int64 > T.high.int64:
            quit "BAD: ev:" & $ev & " c:" & $c

        f.values[se_start] += T(c)
        #echo f.values.len, " ", se_stop, " ", ev
        f.values[se_stop] -= T(c)

  if result and cumsum:
    for f in funs:
      cumulative_sum(f.values)

proc mismatchfun(aln:Record, posns: var seq[mrange]) =
  ## require 4 or more total mismatch "events". NM:i counts a single
  ## insertion of $N bases as $N mismatches. This subtracts so that each
  ## insertion or deletion only counts as a single "mismatch"
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  var nm = tag[int](aln, "NM")
  if nm.isNone or nm.get < 5: return
  var nmi = nm.get
  for c in aln.cigar:
    if c.op == CigarOp.insert or c.op == CigarOp.deletion:
      nmi -= max(1, c.len - 1)
      if nmi < 5: return
  posns.add((aln.start, aln.stop, 1))

proc softfun*(aln:Record, posns:var seq[mrange]) =
  ## softfun an example of a `fun` that can be sent to `genoiser`.
  ## it sets positions where there are soft-clips
  let f = aln.flag
  if f.unmapped or f.qcfail or f.dup: return
  let cig = aln.cigar
  if cig.len == 1: return
  var pos = aln.start

  for op in cig:
    if op.op == CigarOp.soft_clip or op.op == CigarOp.hard_clip:
      # for this function, we want the exact break-points, not the span of the event,
      # so we increment the position and the one that follows it.
      posns.add((pos, pos+1, 1))
    if op.consumes.reference:
      pos += op.len

const endDist = 4 # don't count bases within this many bases of the end of the read.

proc eventfun*(aln:Record, posns:var seq[mrange]) =
  ## eventfun is an example of a `fun` that can be sent to `genoiser`.
  ## it sets positions where there are soft-clips, hard-clips, insertions, or deletions.
  #if aln.mapping_quality == 0: return
  let f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  let cig = aln.cigar
  if cig.len == 1: return
  var pos = aln.start

  let delta = if f.reverse: -1 else: 1

  for op in cig:
    case op.op:
    of CigarOp.soft_clip, CigarOp.hard_clip:
      posns.add((pos, pos+1, delta))
    of CigarOp.insert:
      # for this function, we want the exact break-points, not the span of the event,
      # so we increment the position and the one that follows it.
      if pos - aln.start > endDist and aln.stop - pos > endDist:
        posns.add((pos, pos+1, delta))
    of CigarOp.deletion:
      if pos - aln.start > endDist and aln.stop - pos > endDist:
        posns.add((pos, pos+1, delta))
        posns.add((pos+op.len, pos+op.len+1, delta))
    else:
      discard
    if op.consumes.reference:
      pos += op.len

proc concordant(aln:Record): bool {.inline.} =
  if aln.tid != aln.mate_tid: return false
  # TODO: make these data-driven
  if aln.isize.abs > 1000: return false
  if aln.isize.abs < 20: return false
  var f = aln.flag
  # check we have +- orientation.
  return f.reverse != f.mate_reverse and aln.start > aln.mate_pos == f.reverse


proc read_line(hf: ptr htsFile, kstr: ptr kstring_t, check_ok: proc(depth:int, value:int):bool, idx:int): crange {.inline, thread.} =
  
  var iv = crange(idx:idx)
  while true:
    if hts_getline(hf, cint(10), kstr) <= 0:
      return crange(idx:idx)

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
  return crange(idx:idx)

type sample_count {.shallow.} = object
  arr:seq[int8]
  chrom:string

proc get_chrom_len(f:string, fai:Fai): tuple[chrom:string, length:int] =
  var fh = hts_open(f.cstring, "r")
  var kstr = kstring_t(l:0, m:0, s:nil)
  proc sample_ok(depth:int, value:int): bool = true

  var v = read_line(fh, kstr.addr, sample_ok, 0)
  discard hts_close(fh)
  return (v.chrom, fai.chrom_len(v.chrom))


proc sample_counter(fs: seq[string], sample_checker: string, chrom_len:int, s: ptr sample_count): bool {.thread.}=
  ## this is the 'worker' that gets run in parallel. It takes a batch of `fs` files. and updates the s[].arr
  ## it can not allocate memory related to `s` or the garbage collector could clean it up.
  ## using batches amortizes the cost of the zeroMem (which is not that high).
  var sc = s[]

  if sc.arr.len != chrom_len + 1:
    stderr.write_line "error: cant set new length"
    quit 2
  else:
    zeroMem(sc.arr[0].addr, sizeof(sc.arr[0]) * sc.arr.len)

  if len(fs) > 127:
    stderr.write_line "error: cant batch more thatn 127 samples because of overflow"
    quit 2

  var ex = expression(sample_checker)
  if ex.error() != 0:
    stderr.write_line "error: with expression"
    quit "error with expression"
  var vc = "value".cstring
  var dc = "depth".cstring

  proc sample_ok(depth:int, value:int): bool =
    discard ke_set_int(ex.ke, vc, value)
    discard ke_set_int(ex.ke, dc, depth)
    result = ex.bool
    if ex.error() != 0:
      quit format("expresion error with value:$# depth:$#", value, depth)
  var kstr = kstring_t(l:0, m:0, s:nil)
  # TODO: see why 48119899 is showing up in output.

  for f in fs:
    var fh = hts_open(f.cstring, "r")
    var v = read_line(fh, kstr.addr, sample_ok, 0)
    if v.chrom != "" and sc.chrom != v.chrom:
      quit "got different chromosome values:" & sc.chrom & ", " & v.chrom
    if v.chrom == "":
      stderr.write_line "no lines found in:" & f & " (this can happen for small chromosomes)"
    while v.stop != 0:
      when defined(debug):
        if sc.arr.len == 0:
          stderr.write_line "sc nil"

        if v.stop >= sc.arr.len:
          quit "got stop > chromosome length"
      sc.arr[v.start].inc
      sc.arr[v.stop].dec
      v = read_line(fh, kstr.addr, sample_ok, 0)
    discard hts_close(fh)

  s[] = sc
  result = true


iterator aggregator*(fns: seq[string], sample_checker: string, fai:Fai, nthreads:int): mchrom =
  var threads = min(nthreads, len(fns))

  var (chrom, chrom_len) = get_chrom_len(fns[0], fai)
  var responses = newSeq[FlowVarBase](min(threads, fns.len))
  var results = newSeq[sample_count](min(threads, fns.len))
  var accumulated = newSeq[int32](chrom_len + 1)

  var batch_size = 8

  if batch_size * threads > len(fns):
    batch_size = int(len(fns) / threads)
  if batch_size > 127:
    batch_size = 127

  for j in 0..min(responses.len, results.len) - 1:
    results[j].arr = newSeq[int8](chrom_len+1)
    results[j].chrom = chrom
    #stderr.write_line join(fns[j*batch_size..((j + 1) * batch_size - 1)], "\n")
    responses[j] = spawn sample_counter(fns[j*batch_size..((j + 1) * batch_size - 1)], sample_checker, chrom_len, results[j].addr)

  var jobi = responses.len

  while results.len != 0:

    var index = blockUntilAny(responses)
    if index == -1:
      quit "got unexpected value from await"

    var sc = results[index]
    if chrom == "":
      chrom = sc.chrom
    elif chrom != sc.chrom:
        quit "expected all files from same chromosome, got:'" & chrom & "' and '" & sc.chrom & "'"
    for pos, v in sc.arr:
      # we know the values should not overlap.
      when defined(debug):
        if v < -batch_size or v > batch_size: quit "got unexpected value at:" & $pos & " of " & $v
      accumulated[pos] += v.int32

    if (jobi * batch_size) < fns.len:
      var imin = jobi * batch_size
      var imax = min(fns.len - 1, (jobi + 1) * batch_size - 1)
      #stderr.write_line join(fns[imin..imax], "\n")
      responses[index] = spawn sample_counter(fns[imin..imax], sample_checker, chrom_len, results[index].addr)
    else:
      results.del(index)
      responses.del(index)

    jobi += 1

  cumulative_sum(accumulated)
  for m in ranges(accumulated, chrom):
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
  ## depthfun is an example of a `fun` that can be sent to `genoiser`.
  ## it sets reports the depth at each position
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  #refposns(aln, posns)
  posns.add((aln.start, aln.stop, 1))

proc mq0fun*(aln:Record, posns:var seq[mrange]) =
  ## this is an example function that increments all reference locations with mapping-quality 0.
  if aln.mapping_quality != 0: return
  var f = aln.flag
  if f.unmapped or f.qcfail or f.dup: return
  posns.add((aln.start, aln.stop, 1))

proc mqlt60fun*(aln:Record, posns:var seq[mrange]) =
  ## increment anywhere there's a mapping quality less than 60
  if aln.mapping_quality >= 60'u8: return
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  posns.add((aln.start, aln.stop, 1))

proc aggregate_main(argv: seq[string]) =

  let doc = format("""

  Usage: aggregate [options] <per-sample-output>...

Arguments:

  <per-sample-output>      the txt output from genoiser per-sample

Options:
  -f --fasta <path>         indexed fasta reference file required to get chromosome lengths.
  -t --threads <INT>        number of processors to use for parallelization. [default: 12]
  -n --n-samples <INT>      don't output lines where fewer than this many samples have the noise. larger values result in smaller files. [default: 2]
  -e --expr <string>        per-sample filtering expression [default: "(value > 0) & ((value / depth) > 0.1)"]
  -h --help                 show help
  """)

  let args = docopt(doc, argv=argv)
  #GC_disableMarkAndSweep()

  var
    expr:string = ($args["--expr"]).strip(chars={'"'})
    fai:Fai
    n_samples = parseInt($args["--n-samples"])

  if $args["--fasta"] == "nil":
    quit "--fasta is required"
  if not open(fai, $args["--fasta"]):
    quit "invalid --fasta: " & $args["--fasta"]

  var threads = parseInt($args["--threads"])

  var ex = expression(expr)
  if ex.error() != 0:
    quit "error with expression"

  for iv in aggregator(@(args["<per-sample-output>"]), expr, fai, threads):
    if iv.count < n_samples: continue
    stdout.write_line iv.chrom & "\t" & intToStr(iv.start) & "\t" & intToStr(iv.stop) & "\t" & intToStr(iv.count)

proc writefn(a:Fun, depths:Fun, fh:BGZ, chrom:string, start:int, min_depth:int=0, min_value:int=1) =
  for m in mranges(depths.values, a.values):
    if m.depth < min_depth: continue
    if m.value < min_value: continue
    discard fh.write_line chrom & "\t" & intToStr(m.start + start) & "\t" & intToStr(m.stop + start) & "\t" & intToStr(m.depth) & "\t" & intToStr(m.value)

proc myopen(prefix, chrom, event:string): BGZ =

  createDir(prefix & "/" & event)
  # prefix & ".genoiser." & target.name & ".soft.bed.gz")
  var path = prefix / event / chrom & ".genoiser." & event & ".bed.gz"
  open(result, path, "w1")

proc per_sample_main(argv: seq[string]) =

  let doc = format("""

  Usage: per-sample [options] <PREFIX> <BAM-or-CRAM>

Arguments:

  <PREFIX>       output prefix (must be unique for each sample)
  <BAM-or-CRAM>  the alignment file.

Options:

  -f --fasta <fasta>              fasta file for use with CRAM files.
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
  #discard b.set_option(FormatOption.CRAM_OPT_DECODE_MD, 0)

  if b == nil:
    quit "coudn't open bam: " & fbam
  if b.idx == nil:
    quit "coudn't open bam index for " & fbam

  var L = 5000000
  var depths = Fun[int16](values:new_seq[int16](L+1), f:depthfun)
  var softs = Fun[int16](values:new_seq[int16](L+1), f:softfun)
  #var interc = Fun(values:new_seq[int32](L+1), f:interchromosomal)
  #var weirds = Fun(values:new_seq[int32](L+1), f:weird)
  #var misms = Fun(values:new_seq[int32](L+1), f:mismatchfun)
  #var events = Fun(values:new_seq[int32](L+1), f:eventfun)
  #var mq0 = Fun(values:new_seq[int32](L+1), f:mq0fun)
  #var mqlt60 = Fun(values:new_seq[int32](L+1), f:mqlt60fun)
  #var noise = Fun[int16](values:new_seq[int16](L+1), f:noisefun)
  #var fns = @[depths, softs, weirds, misms, events, interc, mq0, mqlt60]
  var fns = @[depths, softs]

  for target in b.hdr.targets:

    var fhs: seq[BGZ]

    #stderr.write_line target.name
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

      if genoiser[int16](b, fns, target.name, start, stop):
        if fhs.len == 0:
          fhs = @[
            # NOTE: fragile. make sure these are same orders as fns array above.
            myopen(prefix, target.name, "soft"),
            #myopen(prefix, target.name, "weird"),
            #myopen(prefix, target.name, "mismatches"),
            #myopen(prefix, target.name, "event"),
            #myopen(prefix, target.name, "interc"),
            #myopen(prefix, target.name, "mq0"),
            #myopen(prefix, target.name, "mqlt60"),
            #myopen(prefix, target.name, "noise"),
          ]

        #writefn(softs, depths, fhs[0], target.name, start, min_depth=min_depth, min_value=min_value)
        #writefn(weirds, depths, fhs[1], target.name, start, min_depth=min_depth, min_value=min_value)
        #writefn(misms, depths, fhs[2], target.name, start, min_depth=min_depth, min_value=6)
        #writefn(events, depths, fhs[3], target.name, start, min_depth=min_depth, min_value=min_value)
        #writefn(interc, depths, fhs[4], target.name, start, min_depth=min_depth, min_value=min_value)
        #writefn(mq0, depths, fhs[5], target.name, start, min_depth=min_depth, min_value=min_value + 2)
        #writefn(mqlt60, depths, fhs[6], target.name, start, min_depth=min_depth, min_value=min_value + 2)
        writefn(softs, depths, fhs[0], target.name, start, min_depth=min_depth, min_value=min_value)
        #writefn(noise, depths, fhs[1], target.name, start, min_depth=min_depth, min_value=min_value + 2)

        for f in fns:
          zeroMem(f.values[0].addr.pointer, f.values.len * sizeof(f.values[0]))

      start += L

    for f in fhs:
      discard f.close()


var progs = newTable[string, proc (argv: seq[string])]()
progs["per-sample"] = per_sample_main
progs["aggregate"] = aggregate_main

var helps = newTable[string, string]()
helps["per-sample"] = "calculate various metrics for a per-sample BAM/CRAM"
helps["aggregate"] = "aggregate many sample outputs of per-sample into single file"
#var helps = {
#   "per-sample": "calculate various metrics for a per-sample BAM/CRAM",
#   "aggregate": "aggregate many sample outputs of per-sample into single file"
#}.toTable

const version = "0.2.6"

proc main() =

  var args = commandLineParams()
  if len(args) < 1 or (not progs.contains(args[0])):
    var hkeys = toSeq(keys(helps))
    sort(hkeys, proc(a, b: string): int =
      if a < b: return -1
      else: return 1
      )
    echo format("\ngenoiser utility programs.\nversion: $#\n", version)

    for k in hkeys:
      echo format("	• $1: $2", k & repeat(" ", 20 - len(k)), helps[k])
    echo ""
  else:
    var p = args[0]; args.delete(0)
    progs[p](args)

when isMainModule:
  main()
