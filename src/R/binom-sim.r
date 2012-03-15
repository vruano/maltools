#!/usr/bin/env Rscript

KLDist = function(p,q) {
  po = sort(p);
  qo = sort(q);
  ao = sort(unique(c(p,q)));
  i = 1;
  j = 1;
  lp = length(po);
  lq = length(qo);
  pqDist = 0;
  qpDist = 0;
  for (k in 1:length(ao)) {
    aon = ao[k];
    while (i <= lp && po[i] < aon) {
      i = i + 1;
    }
    while (j <= lq && qo[j] < aon) {
      j = j + 1;
    }
    qfreq = j;
    while (j <= lq && qo[j] == aon) {
      j = j + 1;
    }
    qfreq = j - qfreq;
    pfreq = i
    while (i <= lp && po[i] == aon) {
       i = i + 1;
    }
    pfreq = i - pfreq;
    if (qfreq == 0) qfreq = 0.5;
    if (pfreq == 0) pfreq = 0.5;
    pqDist = pqDist  + (pfreq/lp) * log((pfreq/lp)/(qfreq/lq));
    qpDist = qpDist  + (qfreq/lq) * log((qfreq/lq)/(pfreq/lp));
  }
  c(pqDist,qpDist);
}


in.filename = commandArgs()[6];

reps = 100;
lines = 10;
in.file = file(in.filename);
in.lines = readLines(in.file,n=lines);
close(in.file);


processLine = function(l) {
  r = as.numeric(strsplit(l,",")[[1]]);
  avg.rate = mean(10 ^ (-r/10))
  avg.sim = rbinom(reps,length(r),avg.rate);
  exact.sim = sapply(c(1:reps),function (rr) {
    sim.values = sapply(r,function(x) { rbinom(1,1,10 ^ (-x/10)) });
    sum(sim.values);
  });
  ts = ks.test(avg.sim,exact.sim);
  kl = KLDist(avg.sim,exact.sim);
  c(-10 * log(avg.rate,10),ts$p.value,kl,length(r),quantile(avg.sim,prob=0.99),quantile(exact.sim,prob=0.99));
};

in.means = lapply(in.lines,function(x) { processLine (x) });

print(in.means);

