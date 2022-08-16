# reference: https://github.com/ryantibs/best-subset/
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(".")

library(ggplot2)
library(dplyr)

result.file <- "results/tests07.csv"
fig.dir <- "results/fig07"
if (!dir.exists(fig.dir)) dir.create(fig.dir)

dat <- read.csv(result.file, header = T)
dat <- within(dat, {
  beta.type <- as.factor(beta.type)
  rho <- as.factor(rho)
  levels(beta.type) = paste("Beta-type", levels(beta.type))
  levels(rho) = paste("Correlation", levels(rho))
})

n.candidates <- c(100, 200, 500, 1000)
p.candidates <- c(10, 100, 1000, 10000, 30000)
rho.candidates <- c(0, 0.35, 0.7)
s.candidates <- c(5, 10, 20, 50)
beta.type.candidates <- 1:3
snr.candidates <- exp(seq(log(0.05),log(6),length=10))

for (n0 in n.candidates) {
  for (p0 in p.candidates) {
    
    file.name <- paste("n", n0, "p", p0, sep = "")
    pdf(sprintf("%s/%s.pdf", fig.dir, file.name), onefile = TRUE)
    
    for (s0 in s.candidates) {
      for (what in c("obj", "runtime")) {
        df <- dat %>% 
          filter(n == n0, p == p0, s == s0) %>% 
          rename(x = snr, y = what) %>% 
          group_by(n, p, rho, s, beta.type, x, method) %>% 
          summarise(y.avg = mean(y, na.rm=TRUE), y.se = sd(y, na.rm = TRUE))
        if (nrow(df) == 0) next
        xlab <- "Signal-to-noise Ratio"
        ylab <- switch(what,
                       obj="Objective (NLL + penalty)",
                       runtime="Elapsed Time")
        # ylim <- range(df$y.avg-df$y.se, df$y.avg+df$y.se)
        main <- paste0("n=",n0,", p=",p0,", s=",s0, ": \t", what)
        
        gp <- ggplot(df, aes(x=x, y=y.avg, color=method)) + 
          xlab(xlab) + ylab(ylab) + 
          geom_line(lwd=1) + geom_point(pch=19) + 
          theme_bw() + theme(legend.position = "bottom") + 
          ggtitle(main) + 
          geom_errorbar(aes(ymin = y.avg - y.se, ymax = y.avg + y.se), width=0.02) + 
          facet_grid(rows=vars(beta.type), cols=vars(rho))
      
        snr.breaks <- round(exp(seq(from=min(log(dat$snr)),
                                    to=max(log(dat$snr)),length=4)),2)
        gp = gp + scale_x_continuous(trans="log", breaks=snr.breaks)
        
        print(gp)
      }
    }
    dev.off()
  }
}





#============= For Debug Purposes ===============#
# gps <- vector(mode = "list", length = length(beta.type.candidates) * length(rho.candidates))
# i <- 1
# for (beta.type in beta.type.candidates) {
#   for (rho in rho.candidates) {
#     d <- df %>% filter(beta.type==beta.type, rho==rho)
#     ylim <- range(df$y.avg-df$y.se, df$y.avg+df$y.se)
#     gp <- ggplot(df, aes(x=x, y=y.avg, color=method)) + 
#       xlab(xlab) + ylab(ylab) + coord_cartesian(ylim=ylim) +
#       geom_line(lwd=1) + geom_point(pch=19) + 
#       theme_bw() + theme(legend.position = "bottom") + 
#       ggtitle(main) + 
#       geom_errorbar(aes(ymin = y.avg - y.se, ymax = y.avg + y.se), width=0.02)
#     gps[[i]] <- gp
#   }
# }
# 
# library(gridExtra)
# n <- length(gps)
# do.call("grid.arrange", c(gps, ncol=3, nrow=4))


# df.dummy <- df %>% group_by(rho, beta.type) %>%
#   summarise(y.min = min(y.avg-y.se, na.rm=T), y.max = max(y.avg+y.se, na.rm=T), 
#             x=NA)
# df.dummy.min <- df.dummy %>% select(rho, beta.type, y=y.min, x)
# df.dummy.max <- df.dummy %>% select(rho, beta.type, y=y.max, x)
# df.dummy <- bind_rows(df.dummy.min, df.dummy.max)
# df.dummy <- merge(as.data.frame(df.dummy), 
#                   as.data.frame(as.data.frame(df)%>% select(rho, beta.type, method) %>% distinct(rho, beta.type, method)))
# # df.dummy
# 
# gp <- ggplot(df, aes(x=x, y=y.avg, color=method)) + 
#   xlab(xlab) + ylab(ylab) + 
#   geom_blank(data=df.dummy, aes(x=x, y=y, color=method)) +
#   facet_grid(rows=vars(beta.type), cols=vars(rho), scales = "free") + 
#   geom_line(lwd=1) + geom_point(pch=19) + 
#   theme_bw() + theme(legend.position = "bottom") + 
#   ggtitle(main) + 
#   geom_errorbar(aes(ymin = y.avg - y.se, ymax = y.avg + y.se), width=0.02)
