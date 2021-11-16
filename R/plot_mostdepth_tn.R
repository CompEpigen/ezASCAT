#'Plot results from mosdepth output for Tumor/Normal pairs
#' @param t_bed mosdepth output from tumor
#' @param n_bed mosdepth output from normal
#' @param plot_file basename for output plot
#' @param col Colors. Default c("#95a5a6", "#7f8c8d")
#' @export

plot_mosdepth_tn = function(t_bed = NULL, n_bed = NULL, plot_file = NULL, col = NULL){
  
  col = c("#95a5a6", "#7f8c8d")
  
  contigs = c(1:22, "X", "Y", paste0("chr", 1:22), "chrX", "chrY")
  
  dat = lapply(X = c(t_bed, n_bed), function(x){
    x = data.table::fread(input = x)
    colnames(x) = c("chr", "start", "end", "doc")
    x[, chr := gsub(pattern = "chr", replacement = "", x = chr)]
    colnames(x)[1:3] = c("Chromosome", "Start_Position", "End_Position")
    x = x[Chromosome %in% contigs]
    x
  })
  
  
  names(dat) = c("tumor", "normal")
  dat = merge(dat$tumor, dat$normal, by = c("Chromosome", "Start_Position", 'End_Position'), suffixes = c("_t", "_n"))
  
  dat_xy = dat[Chromosome %in% c('X', 'Y', 'chrX', 'chrY')]
  datn = dat[!Chromosome %in% c('X', 'Y', 'chrX', 'chrY')]
  datn[, Chromosome := as.numeric(as.character(Chromosome))]
  datn = datn[order(Chromosome, Start_Position)]
  dat = rbind(datn, dat_xy)
  
  #Get chr lengths
  chr.lens.dt = dat[,max(End_Position, na.rm = TRUE), .(Chromosome)]
  chr.lens = chr.lens.dt$V1
  names(chr.lens) = chr.lens.dt$Chromosome
  
  map_ratio = sum(dat$doc_t, na.rm = TRUE)/sum(dat$doc_n, na.rm = TRUE)
  message("Coverage ration T/N: ", round(map_ratio, digits = 3))
  dat$doc_n = dat$doc_n * map_ratio
  
  dat[, logR := log2(doc_t+1) - log2(doc_n+1)]
  
  cols = rep(x = c("#95a5a6", "#7f8c8d"), length(chr.lens))
  
  
  all_depth_spl = split(dat, dat$Chromosome)
  all_depth_spl = all_depth_spl[names(chr.lens)]
  
  seg.spl.transformed = all_depth_spl[[1]]
  if (nrow(seg.spl.transformed) > 0) {
    seg.spl.transformed$Start_Position_updated = seg.spl.transformed$Start_Position
    seg.spl.transformed$End_Position_updated = seg.spl.transformed$End_Position
  }
  chr.lens.sumsum = cumsum(as.numeric(chr.lens))
  for (i in 2:length(all_depth_spl)) {
    x.seg = all_depth_spl[[i]]
    if (nrow(x.seg) > 0) {
      x.seg$Start_Position_updated = x.seg$Start_Position + 
        chr.lens.sumsum[i - 1]
      x.seg$End_Position_updated = x.seg$End_Position + 
        chr.lens.sumsum[i - 1]
    }
    seg.spl.transformed = rbind(seg.spl.transformed, x.seg, 
                                fill = TRUE)
  }
  
  all_depth_spl = split(seg.spl.transformed, seg.spl.transformed$Chromosome)
  all_depth_spl = all_depth_spl[names(chr.lens)]
  
  rm(seg.spl.transformed)
  
  cols = rep(x = c("#95a5a6", "#7f8c8d"), length(chr.lens))
  
  png(filename = paste0("xx", ".png"), width = 1024, height = 600, bg = "white")
  plot(NA, xlim = c(0, sum(chr.lens)), ylim = c(-3, 3), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  temp = lapply(seq_along(all_depth_spl), function(idx){
    x = all_depth_spl[[idx]]
    points(x$Start_Position_updated, x$logR, pch = "-", col = cols[idx])
    # rect(xleft = x[, Start_Position_updated][1], ybottom = log2(med_cov),
    #      xright = x[,End_Position_updated][nrow(x)], ytop = log2(med_cov))
  })
  abline(v = cumsum(as.numeric(chr.lens)), lty = 2)
  axis(side = 1, at = cumsum(as.numeric(chr.lens)), labels = names(chr.lens))
  axis(side = 2, at = seq(-3, 3, 1), las = 2)
  title(main = "DOC Median centered")
  dev.off()
}