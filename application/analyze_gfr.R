par.old <- par(no.readonly = TRUE)

load(file.path("application", "result_gfr.Rdata"))
load(file.path("application", "result_gfr_olcov.Rdata"))

fig_path <- "application"

tgrid1 <- seq(1,7,length.out=101)
tgrid2 <- seq(1,7,length.out=50)

check_points <- c(20, 100, 500)
q <- 3
nBlock <- 200
nums = check_points * nBlock
num_labels = c()
for (idx in seq_along(nums)) {
  if (nums[idx] < 1000) {
    num_labels[idx] = as.character(nums[idx])
  } else {
    num_labels[idx] = paste0(as.character(nums[idx] %/% 1000), ",000")
  }
}

# flip the FPC if most values are negative
flip_idx = colMeans(PhiAvgAll[,,check_points[3]]) < 0
PhiAvgAll[,flip_idx,] = -PhiAvgAll[,flip_idx,]
PhiEstAll.ll[,flip_idx,] = -PhiEstAll.ll[,flip_idx,]

phi_lims = sapply(1:q, \(k) {
  range(c(PhiAvgAll[-(1:4),k,]) |> na.omit())
})

fve = lambda.avg / sum(lambda.avg)

saveRDS(PhiAvgEval, file.path("application","eigf_gfr.rds"))
saveRDS(lambda.avg, file.path("application","eigval_gfr.rds"))

setEPS()
postscript(file=file.path(fig_path, "gfr_eigfun.eps"), width=7, height=6.2)
npanel = prod(c(q, length(check_points)))
layout(mat = matrix(c(1:npanel, rep(npanel+1, q)), byrow=TRUE,
                    nrow=length(check_points)+1, ncol=q),
       heights = c(0.3,0.3,0.3,0.1))
par(mar=c(2,4.5,2,1.5), oma=c(2,0,1,0))
for (idx in seq_along(check_points)) {
  for (k in 1:q) {
    plot(tgrid1, PhiAvgAll[,k,check_points[idx]+1], type="l", col="#A21B2D", lty=1, lwd=3,
         xlab="", ylab="", main=NULL, ylim = phi_lims[,k])
    lines(tgrid2, PhiEstAll.ll[,k,check_points[idx]], lty=4, col="#1A73A0", lwd=3)
    if (k==1) {
      mtext(paste0("n=",num_labels[idx]), side = 2, line = 3)
    }
    if (idx==1) {
      # mtext(paste(numtext, "FPC,", "FVE=", round(fve[k]*100,1), "%"), side = 3, line = 1)
      mtext(paste0("FPC ", k, ", ", round(fve[k]*100,1), "%"), side = 3, line = 1)
    }
    if (idx==3 && k == 2) {
      mtext("Years After Transplant", side = 1, line = 3)
    }
  }
}
par(mar=c(0,2,1,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom", legend=c("OnlineFPCA-RSGD", "OnlineCov"), inset=c(0,0),
       lty=c(1,4), col=c("#A21B2D", "#1A73A0"), lwd=2,
       cex=1.2, horiz=TRUE, bty="n")
par(par.old)
dev.off()

setEPS()
postscript(file=file.path(fig_path,"gfr_taupath.eps"), width=9, height=5)
with(tau.select, plot.tau_path2(tau.history, tau.selectId))
dev.off()

setEPS()
postscript(file=file.path(fig_path,"gfr_time.eps"), width = 6, height=5)
matplot((1:718) * 200 / 1000,
  cbind(cumsum(as.numeric(sgd_time[-1])), cumsum(as.numeric(loclin_time))) / 60,
  type="l", xlab="Number of Recipients (Thousand)", ylab="CPU Time (min)",
  col = c("#A21B2D", "#1A73A0"), lwd=3)
legend("topleft", legend=c("OnlineFPCA-RSGD", "OnlineCov"),
       lty=1:2, col = c("#A21B2D", "#1A73A0"), lwd=2)
dev.off()
