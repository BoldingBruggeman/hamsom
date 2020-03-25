c-----------------------------------------------------------------------
c     x,y e.g. phi lam coordinates in grad
c-----------------------------------------------------------------------
      common /xy_coor/  xt(m), yt(n),cosxt(m)	 ! coordinates of T-points
      common /xy_cooru/ xv(m), yu(n)           ! coordinates of uv-points
      common /part_dens/ part_d(m,n,ilo)
      common /begin_xy/ phigrad,yamgrad

