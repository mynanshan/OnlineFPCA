library(fda)

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/onlineFPCA.R")

par.old <- par(no.readonly = TRUE)

load(file.path("application", "result_gfr.Rdata"))
load(file.path("application", "result_gfr_olcov.Rdata"))

fig_path <- "application"

tgrid <- seq(1,7,length.out=31)
tgrid1 <- seq(1,7,length.out=101)
tgrid2 <- seq(1,7,length.out=50)
basis <- create.bspline.basis(c(1, 7), nbasis = 7, norder = 4)

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

# order FPCs by FVEs
ord <- order(lambda.avg, decreasing = TRUE)
PhiAvgEval <- PhiAvgEval[,ord]
PhiAvgAll <- PhiAvgAll[,ord,]
resCI = resCI[,ord,,]
PhiEst.ll <- PhiEst.ll[,ord]
PhiEstAll.ll <- PhiEstAll.ll[,ord,]
lambda.avg <- lambda.avg[ord]
lambdaEst.ll <- lambdaEst.ll[ord]

# flip the FPC if most values are negative
flip_idx = colMeans(PhiAvgAll[,,check_points[3]]) < 0
PhiAvgAll[,flip_idx,] = -PhiAvgAll[,flip_idx,]
flip_idx = colMeans(PhiEstAll.ll[,,check_points[3]]) < 0
PhiEstAll.ll[,flip_idx,] = -PhiEstAll.ll[,flip_idx,]

phi_lims = sapply(1:q, \(k) {
  range(c(PhiAvgAll[-(1:4),k,-(1:check_points[1])]) |> na.omit())
})
all_lims = sapply(1:q, \(k) {
  range(c(PhiAvgAll[-(1:4),k,-(1:check_points[1])],
  PhiEstAll.ll[-(1:4),k,-(1:check_points[1])], resCI[-(1:4),k,,-(1:check_points[1])]) |> na.omit())
})
all_lims <- c(0.5, 1.3, -2.2, 4, -2.2, 4) |> 
  matrix(nrow = 2)

fve = lambda.avg / sum(lambda.avg)

saveRDS(PhiAvgEval, file.path("application","eigf_gfr.rds"))
saveRDS(lambda.avg, file.path("application","eigval_gfr.rds"))

pdf(file=file.path(fig_path, "gfr_eigfun.pdf"), width=7, height=6.2)
npanel = prod(c(q, length(check_points)))
add_ci = TRUE
layout(mat = matrix(c(1:npanel, rep(npanel+1, q)), byrow=TRUE,
                    nrow=length(check_points)+1, ncol=q),
       heights = c(0.3,0.3,0.3,0.1))
par(mar=c(2,4.5,2,1.5), oma=c(2,0,1,0))
for (idx in seq_along(check_points)) {
  for (k in 1:q) {
    phik.est <- PhiAvgAll[,k,check_points[idx]+1]
    phik.olcov <- PhiEstAll.ll[,k,check_points[idx]]
    plot(
      tgrid1, phik.est, type="l", col="#A21B2D", lty=1, lwd=3,
      xlab="", ylab="", main=NULL,
      ylim = if (add_ci) all_lims[,k] else phi_lims[,k]
    )
    lines(tgrid2, phik.olcov, lty=4, col="#1A73A0", lwd=3)
    if (add_ci) {
      phik.l <- resCI[, k, 1, check_points[idx]+1]
      phik.u <- resCI[, k, 3, check_points[idx]+1]
      if (
        inprod(
          smooth.basis(tgrid, resCI[, k, 2, check_points[idx]+1], basis)$fd,
          smooth.basis(tgrid1, phik.est, basis)$fd
        ) < 0
      ) {
        phik.l <- -phik.l
        phik.u <- -phik.u
      }
      polygon(
        c(tgrid, rev(tgrid)),
        c(phik.l, rev(phik.u)),
        col = rgb(0.2, 0.3, 0.9, alpha = 0.2),
        border = NA
      )
    }
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

pdf(file=file.path(fig_path,"gfr_taupath.pdf"), width=7, height=4.5)
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
