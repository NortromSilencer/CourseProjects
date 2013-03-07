!     
! File:   locate.f
! Author: Peisi
!
! Created on March 7, 2013, 12:32 PM
!

!----------------------------------------------------------------------
      Subroutine Locate(XX,N,X,J)
      Implicit Real(8) (A-H,O-Z)
      Dimension XX(N)
      JL=0
      JU=N+1
10    if(JU-JL > 1)then
        JM=(JU+JL)/2
        if((XX(N) >= XX(1)).EQV.(X >= XX(JM)))then
          JL=JM
        else
          JU=JM
        endif
      GOTO 10
      endif
      if(X == XX(1))then
        J=1
      else if(X == XX(N))then
        J=N-1
      else
        J=JL
      endif
      return
      end

