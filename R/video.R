# Copyright (C) 2013 Philipp Benner
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

x        <- 0:100/20
theta    <- (1-tanh(-50:50/20))/2.5+0.1

e        <- new.experiment(c(1.0,1.0),
                           partially.apply(kernel.exponential, 0.8))
sample.x <- 0:10/2
gt       <- new.gt(x, theta)

for (i in 1:200) {
  sample.with.gt(e, sample.x, gt)
  png(filename=sprintf("plot_%03d.png", i),
      width = 1200, height = 800)
  plot(e, x)
  title(main = sprintf("%d samples", i))
  lines(x, theta, 'l', lwd=3, col='red')
  dev.off()
}

# mencoder mf://plot_*.png -mf type=png:fps=4 -ovc lavc -lavcopts
# vcodec=mpeg4 -oac copy -o plot.avi
