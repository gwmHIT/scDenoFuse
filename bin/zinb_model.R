# Original scFusion nine-parameter ZINB; named input and cached numeric vectors.
# No realignment seeds, LSTM or read identities are used. Direction metadata
# passes through to output only; it never enters fitting or P-value simulation.

parse_numbers = function(text) {
  body = substr(text, 2, nchar(text) - 1)
  if (body == "") return(numeric(0))
  as.numeric(strsplit(body, ",", fixed = TRUE)[[1]])
}

cache_vectors = function(data) {
  data$counts = lapply(data$SplitReadCountsByCell, parse_numbers)
  data$expr1 = lapply(data$Gene1ExpressionBySplitCells, parse_numbers)
  data$expr2 = lapply(data$Gene2ExpressionBySplitCells, parse_numbers)
  data
}

# The legacy code serialized transformed expressions with as.character before
# reparsing them in optim. Retain that rounding, not just the mathematical value.
legacy_roundtrip = function(values) {
  parse_numbers(paste0("[", paste(as.character(values), collapse = ", "), "]"))
}

make_zinb_objective = function(data, gcprelist, numcell) {
  # These values do not depend on the nine fitted parameters. Prepare them
  # once, after the unchanged background selection/expression transformation.
  rows = lapply(seq_len(nrow(data)), function(i) {
    cc = data$counts[[i]]
    expr1 = data$expr1[[i]]
    expr2 = data$expr2[[i]]
    if (sum(is.na(expr1)) > 0 | sum(is.na(expr2)) > 0 | sum(is.na(cc)) > 0) return(NULL)
    gc = (data$GC1[i] + data$GC2[i]) / 2
    list(cc = cc, expr1 = expr1, expr2 = expr2,
         aveexpr1 = mean(expr1), aveexpr2 = mean(expr2),
         gcpre = log(gcprelist[round(gc * 1000) + 1]))
  })
  function(par) {
    loglik = 0
    lambda = exp(par[5])
    coefficient1 = exp(par[3])
    coefficient2 = exp(par[4])
    for (row in rows) {
      if (is.null(row)) next
      cc = row$cc
      expr1 = row$expr1
      expr2 = row$expr2
      gcpre = row$gcpre
      linpart = par[6] + par[7] * gcpre + par[8] * row$aveexpr1 + par[9] * row$aveexpr2
      zeroprob = min(exp(linpart) / (1 + exp(linpart)), 0.9999999999999)
      # Keep the original scalar addition/subtraction order, not a vector sum.
      for (j in seq_along(cc)) {
        mu = par[1] + par[2] * gcpre + coefficient1 * expr1[j] + coefficient2 * expr2[j]
        mu = exp(mu)
        loglik = loglik + dnbinom(as.numeric(cc[j]), 1 / lambda,
                    prob = 1 / (1 + mu * lambda), log = TRUE) + log(1 - zeroprob)
        loglik = loglik - log(1 - dnbinom(0, 1 / lambda, prob = 1 / (1 + mu * lambda)))
      }
      loglik = loglik + (numcell - length(cc)) * log(zeroprob)
      loglik = loglik - log(1 - zeroprob ^ numcell - numcell * (1 - zeroprob) * zeroprob ^ (numcell - 1))
      if (is.nan(loglik)) break
    }
    -loglik
  }
}

calc_zinb = function(par, data, gcprelist, numcell) {
  make_zinb_objective(data, gcprelist, numcell)(par)
}

calc_pvalues = function(par, data, gcprelist, numcell, with_log = FALSE) {
  pv = rep(0, nrow(data))
  logpv = rep(0, nrow(data))
  for (i in seq_len(nrow(data))) {
    cc = data$counts[[i]]
    if (sum(cc) <= 3 | data$CellSupport[i] <= 1) {
      pv[i] = 1
      next
    }
    expr1 = data$expr1[[i]]
    expr2 = data$expr2[[i]]
    if (sum(is.na(expr1)) > 0 | sum(is.na(expr2)) > 0 | sum(is.na(cc)) > 0) {
      pv[i] = 1
      next
    }
    aveexpr1 = mean(expr1) / 2
    aveexpr2 = mean(expr2) / 2
    gc = (data$GC1[i] + data$GC2[i]) / 2
    gcpre = log(gcprelist[round(gc * 1000) + 1])
    linpart = par[6] + par[7] * gcpre + par[8] * aveexpr1 + par[9] * aveexpr2
    zeroprob = min(exp(linpart) / (1 + exp(linpart)), 0.9995)
    for (j in seq_along(cc)) cc[j] = min(15, cc[j])
    totalread = sum(cc)
    mu = par[1] + par[2] * gcpre + exp(par[3]) * 2 * aveexpr1 + exp(par[4]) * 2 * aveexpr2
    truemu = exp(mu)
    mu = max(truemu, 0.1)
    lambda = exp(par[5])
    truevar = truemu + truemu ^ 2 * lambda
    lambda = max(0.001, (truevar - mu) / (mu ^ 2))
    simuset = rnbinom(10000, 1 / lambda, prob = 1 / (1 + mu * lambda))
    simuset = simuset[simuset > 0]
    whilecount = 0
    while (length(simuset) < 2000) {
      newsimuset = rnbinom(10000, 1 / lambda, prob = 1 / (1 + mu * lambda))
      newsimuset = newsimuset[newsimuset > 0]
      simuset = c(simuset, newsimuset)
      whilecount = whilecount + 1
      if (whilecount > 100) simuset = c(simuset, 1)
    }
    simuset[simuset >= 15] = 15
    totalreadsample = c()
    for (j in seq_len(100)) {
      nonzeronum = rbinom(1, numcell, 1 - zeroprob)
      usedreadnum = sample(simuset, nonzeronum, replace = TRUE)
      totalreadsample = c(totalreadsample, sum(usedreadnum))
    }
    pv[i] = pnorm(totalread, mean(totalreadsample), sd(totalreadsample), lower.tail = FALSE)
    if (with_log) {
      logpv[i] = pnorm(totalread, mean(totalreadsample), sd(totalreadsample),
                      lower.tail = FALSE, log.p = TRUE)
    }
  }
  if (with_log) return(list(pvalues = pv, log_pvalues = logpv))
  pv
}

gc_background = function(data) {
  data = data[!(pmax(data$GC1, data$GC2) <= 0.05 | pmin(data$GC1, data$GC2) >= 0.95), ]
  gc = round((data$GC1 + data$GC2) / 2, digits = 3)
  gclist = unique(round(sort(unique(gc)), digits = 3))
  chimexpr = gclist
  for (i in seq_along(gclist)) {
    use = which(gc == gclist[i])
    exprtotal = 0
    count = 0
    for (j in use) {
      cc = data$counts[[j]]
      if (length(cc) >= 1) {
        exprtotal = exprtotal + sum(cc)
        count = count + length(cc)
      }
    }
    chimexpr[i] = exprtotal / count
  }
  lmres = lm(chimexpr ~ splines::bs(gclist, df = 5, intercept = TRUE))
  basepre = predict(lmres, data.frame(gclist = gclist))
  diff = abs(chimexpr - basepre)
  threshold = mean(diff) + 1.5 * sqrt(var(diff))
  for (i in seq_along(gclist)) {
    if (chimexpr[i] - basepre[i] > threshold) chimexpr[i] = basepre[i]
  }
  lmres = lm(chimexpr ~ splines::bs(gclist, df = 5, intercept = TRUE))
  curve = predict(lmres, data.frame(gclist = seq(0, 1, length.out = 1001)))
  curve[curve <= 0] = 0.01
  curve
}

prepare_scoring = function(data) {
  data = data[!(pmax(data$GC1, data$GC2) <= .05 | pmin(data$GC1, data$GC2) >= .95), ]
  keep = rep(TRUE, nrow(data))
  data$totalreadcount = rep(0, nrow(data))
  for (i in seq_len(nrow(data))) {
    cc = data$counts[[i]]
    data$totalreadcount[i] = sum(cc)
    smallindex = which(cc <= 2)
    bigindex = which(cc > 2)
    if (max(cc) < 2 | length(cc) < 2) keep[i] = FALSE
    usesmallcellnum = floor(data$SpanningSupportingCellCount[i] * 1.1 + 2)
    if (sum(cc <= 2) <= usesmallcellnum) {
      usesmallcell = smallindex
    } else {
      usesmallcell = sample(smallindex, usesmallcellnum)
    }
    usedindex = union(bigindex, usesmallcell)
    if (!keep[i]) next
    gene1expr = data$expr1[[i]]
    gene2expr = data$expr2[[i]]
    if (min(min(gene1expr), min(gene2expr)) < 0.0001) {
      keep[i] = FALSE
      next
    }
    gc1 = data$GC1[i]
    gc2 = data$GC2[i]
    data$GC1[i] = min(gc1, gc2)
    data$GC2[i] = max(gc1, gc2)
    data$counts[[i]] = legacy_roundtrip(cc[usedindex])
    data$expr1[[i]] = legacy_roundtrip(log(1 + pmin(gene1expr, gene2expr) / 100)[usedindex])
    data$expr2[[i]] = legacy_roundtrip(log(1 + pmax(gene1expr, gene2expr) / 100)[usedindex])
  }
  data[keep, ]
}

prepare_background = function(data) {
  keep = rep(TRUE, nrow(data))
  lastpos1 = c(-1)
  lastpos2 = c(-1)
  cellquantiles = quantile(data$CellSupport)
  highthreshold = cellquantiles[3] + cellquantiles[4] - cellquantiles[2]
  for (i in seq_len(nrow(data))) {
    cc = data$counts[[i]]
    a = runif(1)
    if (a > mean(cc) ^ .1 * 70000 / nrow(data)) {
      keep[i] = FALSE
      next
    }
    if (max(cc) < 2 && length(cc) <= 5) keep[i] = FALSE
    if (data$CellSupport[i] > highthreshold && nrow(data) < 150) keep[i] = FALSE
    for (k in seq_along(lastpos1)) {
      if ((data$Breakpoint1[i] == lastpos1[k] && data$Breakpoint2[i] == lastpos2[k]) ||
          (data$Breakpoint1[i] == lastpos2[k] && data$Breakpoint2[i] == lastpos1[k])) {
        keep[i] = FALSE
        break
      }
    }
    if (keep[i]) {
      lastpos1 = c(lastpos1, data$Breakpoint1[i])
      lastpos2 = c(lastpos2, data$Breakpoint2[i])
      if (length(lastpos1) > 100) {
        lastpos1 = lastpos1[2:101]
        lastpos2 = lastpos2[2:101]
      }
    }
    if (!keep[i]) next
    # Retain the first-pass legacy in-place GC ordering, including its order
    # of assignment, so this refactor does not silently change model estimates.
    data$GC1[i] = min(data$GC1[i], data$GC2[i])
    data$GC2[i] = max(data$GC1[i], data$GC2[i])
    gene1expr = data$expr1[[i]]
    gene2expr = data$expr2[[i]]
    data$expr1[[i]] = legacy_roundtrip(log(1 + pmin(gene1expr, gene2expr) / 100))
    data$expr2[[i]] = legacy_roundtrip(log(1 + pmax(gene1expr, gene2expr) / 100))
    if (data$Chrom1[i] == data$Chrom2[i] && abs(data$Breakpoint2[i] - data$Breakpoint1[i]) < 200000) keep[i] = FALSE
    if (min(min(gene1expr), min(gene2expr)) < .0001) keep[i] = FALSE
  }
  data = data[keep, ]
  keep = rep(TRUE, nrow(data))
  for (i in seq_len(nrow(data))) {
    a = runif(1)
    keep[i] = !(a > mean(data$counts[[i]]) ^ .1 * 700 / nrow(data))
  }
  data[keep, ]
}

read_named = function(path) {
  cache_vectors(read.delim(path, header = TRUE, quote = "", comment.char = "", check.names = FALSE,
                          stringsAsFactors = FALSE))
}

main = function() {
  options(warn = -1)
  argv = commandArgs(trailingOnly = TRUE)
  args = as.list(argv[seq(2, length(argv), by = 2)])
  names(args) = argv[seq(1, length(argv), by = 2)]
  numcell = as.integer(args[["--cells"]])
  data = read_named(args[["--candidates"]])
  background = read_named(args[["--background"]])
  set.seed(1122)
  gcprelist = gc_background(data)
  data = prepare_scoring(data)
  background = prepare_background(background)
  if (nrow(background) == 0) stop("No eligible background candidates for ZINB fitting")
  if (is.null(args[["--reuse-parameters"]])) {
    message("Estimating original ZINB parameters (sequential, seed 1122)")
    objective = make_zinb_objective(background, gcprelist, numcell)
    optres = optim(rnorm(9, mean = -.4, sd = .1), objective,
                   control = list(maxit = 10000))
    rm(objective)
  } else {
    saved = new.env()
    load(args[["--reuse-parameters"]], envir = saved)
    if (!identical(numcell, saved$numcell) || !identical(gcprelist, saved$gcprelist)) {
      stop("Saved model cells or GC curve differ from the current input")
    }
    # The original Nelder-Mead objective is deterministic. Replay preparation
    # and the nine initial normal draws, without repeating optimization.
    invisible(rnorm(9, mean = -.4, sd = .1))
    if (exists("scoring_rng", envir = saved, inherits = FALSE) &&
        !identical(.Random.seed, saved$scoring_rng)) {
      stop("Saved scoring RNG differs from replayed preparation")
    }
    optres = saved$optres
    message("Reusing saved fit; original sequential scoring RNG replayed")
  }
  scoring_rng = .Random.seed
  save(optres, gcprelist, numcell, scoring_rng, file = args[["--parameters"]])
  rm(background)
  message("Calculating original sequential P-values")
  with_log = !is.null(args[["--log-output"]])
  scores = calc_pvalues(optres$par, data, gcprelist, numcell, with_log = with_log)
  pv = if (with_log) scores$pvalues else scores
  spanning = vapply(seq_len(nrow(data)), function(i) {
    # Preserve the original first-pass V6 total, including spanning-only cells.
    # The threshold's R feature separately uses the named per-cell table.
    sum(parse_numbers(data$SpanningReadCountsByCell[i]))
  }, numeric(1))
  output = data.frame(FusionName = paste(data$Gene1, data$Gene2, sep = "--"),
                      CellSupport = data$CellSupport, JunctionReadCount = data$totalreadcount,
                      SpanningFragCount = spanning,
                      Position1 = paste(data$Chrom1, data$Breakpoint1, sep = ":"),
                      Position2 = paste(data$Chrom2, data$Breakpoint2, sep = ":"),
                      Pv23 = pv, CandidateID = data$CandidateID)
  for (field in c("direction1", "direction2", "strand1", "strand2")) {
    output[[field]] = if (field %in% names(data)) data[[field]] else rep("?", nrow(data))
  }
  write.table(output, args[["--output"]], sep = "\t", quote = FALSE, row.names = FALSE)
  if (with_log) {
    log_output = data.frame(CandidateID = data$CandidateID, Pv23 = pv,
                            LogPv23 = scores$log_pvalues)
    write.table(log_output, args[["--log-output"]], sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

if (sys.nframe() == 0) main()
