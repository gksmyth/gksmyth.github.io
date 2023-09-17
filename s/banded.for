       subroutine choleb(n,m,A)
c      LDL' decomposition for banded p.d. matrix
c      input:
c      n - order of matrix
c      m - number of bands, including diagonal
c      A - matrix, for i.ge.j ij element stored in A(i,i-j+1)
c      output:
c      A - D stored in first column, L in columns 2 to m
c
c      Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
c      Last revised 1983.
c
       integer n,m
       double precision A(n,m)
       integer i,j,k
       double precision r,s
c
       do 1 i=1,n
         do 2 j=max0(1,i-m+1),i-1
           s=0.d0
           do 3 k=max0(1,i-m+1),j-1
    3        s=s+A(i,i-k+1)*A(j,j-k+1)
    2      A(i,i-j+1)=A(i,i-j+1)-s
         s=0.d0
         do 4 j=max0(1,i-m+1),i-1
           r=A(i,i-j+1)/A(j,1)
           s=s+r*A(i,i-j+1)
    4      A(i,i-j+1)=r
    1    A(i,1)=A(i,1)-s
c
       return
       end

       subroutine solblt(n,m,L,b)
c      overwrite b with solution of Lx=b
c      L - lower triangular, unit diagonal, order n
c          m-1 bands below the diagonal, stored in
c          nxm matrix with i,jth element in L(i,i-j+1)
c
c      Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
c      Last revised 1983.
c
       double precision L(n,1), b(1),s
       integer n,m,im1,i,j,i1
c
       im1=1
       do 1 i=2,n
          i1=i+1
          if(i.gt.m) im1=im1+1
          s=0.d0
          do 2 j=im1,i-1
             ij1=i1-j
    2        s=s+L(i,ij1)*b(j)
    1     b(i)=b(i)-s
       return
       end

       subroutine solbut(n,m,U,b)
c      overwrite b with solution of Ux=b
c      U - upper triangular, unit diagonal, order n,
c          m-1 bands above the diagonal, stored in
c          nxm array with i,jth element in U(j,j-i+1)
c
c      Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
c      Last revised 1983.
c
       double precision U(n,1),b(1),s
       integer m,n,i,j,i1,im1
c
       im1=n
       do 1 i=n-1,1,-1
          i1=i-1
          if((n-i1).gt.m) im1=im1-1
          s=0.d0
          do 2 j=im1,i+1,-1
    2        s=s+U(j,j-i1)*b(j)
    1     b(i)=b(i)-s
       return
       end
