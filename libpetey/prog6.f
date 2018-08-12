
      program PROG6
c
C  DEMONSTRATE THE USE OF THE SUBROUTINE   LDP  FOR LEAST   
C  DISTANCE PROGRAMMING BY SOLVING THE CONSTRAINED LINE DATA FITTING 
C  PROBLEM OF CHAPTER 23.
C   
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 15, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------   
      integer MX
      parameter(MX = 2)
      integer I, INDEX(3), J, L, MDE, MDGH, ME, MG, MODE, N, NP1
      double precision E(4,MX), F(4), G(3,MX), G2(3,MX)
      double precision H(3), H2(3), ONE, RES, S(MX), SM, T(4)
      double precision W(4), WLDP(21), WORK(2*MX), X(MX)
      double precision Z(MX), ZERO, ZNORM
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
      data T / 0.25d0, 0.5d0, 0.5d0, 0.8d0 /
      data W / 0.50d0, 0.6d0, 0.7d0, 1.2d0 /
C     ------------------------------------------------------------------   
  110 format (' PROG6..  EXAMPLE OF CONSTRAINED CURVE FITTING'/
     * 10x,'USING THE SUBROUTINE LDP.'//
     * 10x,'RELATED INTERMEDIATE QUANTITIES ARE PRINTED.')  
  120 format (/' V =      ',2F10.5/(10X,2F10.5))
  130 format (/' F TILDA =',4F10.5/' S =      ',2F10.5)    
  140 format (/' G TILDA =',2F10.5/(10X,2F10.5))
  150 format (/' H TILDA =',3F10.5) 
  160 format (/' Z =      ',2F10.5) 
  170 format (/' THE COEFICIENTS OF THE FITTED LINE F(T)=X(1)*T+X(2)'/
     * ' ARE X(1) = ',F10.5,' AND X(2) = ',F10.5)  
  180 format (/' THE CONSECUTIVE RESIDUALS ARE'/1X,4(I4,F10.5))
  190 format (/' THE RESIDUALS NORM IS ',F10.5) 
  200 format (/' MODE (FROM LDP) = ',I3,',  ZNORM = ',F10.5)  
C     ------------------------------------------------------------------   
      write (*,110)     
      MDE=4 
      MDGH=3
C   
      ME=4  
      MG=3  
      N=2   
C                DEFINE THE LEAST SQUARES AND CONSTRAINT MATRICES.  
C   
          DO 10 I=1,ME  
          E(I,1)=T(I)   
          E(I,2)= ONE     
   10     F(I)=W(I)     
C   
      G(1,1)= ONE 
      G(1,2)= ZERO 
      G(2,1)= ZERO 
      G(2,2)= ONE 
      G(3,1)=-ONE
      G(3,2)=-ONE
C   
      H(1)= ZERO   
      H(2)= ZERO   
      H(3)=-ONE  
C   
C     COMPUTE THE SINGULAR VALUE DECOMPOSITION OF THE MATRIX, E.
C   
      CALL SVDRS (E, MDE, ME, N, F, 1, 1, S, WORK)
C   
      write (*,120) ((E(I,J),J=1,N),I=1,N)  
      write (*,130) F,(S(J),J=1,N)  
C   
C     GENERALLY RANK DETERMINATION AND LEVENBERG-MARQUARDT  
C     STABILIZATION COULD BE INSERTED HERE.     
C   
C        DEFINE THE CONSTRAINT MATRIX FOR THE Z COORDINATE SYSTEM.  
          DO 30 I=1,MG  
              DO 30 J=1,N   
              SM= ZERO     
                  DO 20 L=1,N   
   20             SM=SM+G(I,L)*E(L,J)   
   30         G2(I,J)=SM/S(J)   
C         DEFINE CONSTRAINT RT SIDE FOR THE Z COORDINATE SYSTEM.
          DO 50 I=1,MG  
          SM= ZERO 
              DO 40 J=1,N   
   40         SM=SM+G2(I,J)*F(J)    
   50     H2(I)=H(I)-SM 
C   
      write (*,140) ((G2(I,J),J=1,N),I=1,MG)    
      write (*,150) H2  
C   
C                        SOLVE THE CONSTRAINED PROBLEM IN Z-COORDINATES.
C   
      CALL LDP (G2,MDGH,MG,N,H2,Z,ZNORM,WLDP,INDEX,MODE)    
C   
      write (*,200) MODE,ZNORM  
      write (*,160) Z   
C   
C                    TRANSFORM BACK FROM Z-COORDINATES TO X-COORDINATES.
          DO 60 J=1,N   
   60     Z(J)=(Z(J)+F(J))/S(J)     
          DO 80 I=1,N   
          SM= ZERO 
              DO 70 J=1,N   
   70         SM=SM+E(I,J)*Z(J)     
   80     X(I)=SM   
      RES=ZNORM**2  
      NP1=N+1   
          DO 90 I=NP1,ME
   90     RES=RES+F(I)**2   
      RES=sqrt(RES)     
C                                    COMPUTE THE RESIDUALS. 
          DO 100 I=1,ME 
  100     F(I)=W(I)-X(1)*T(I)-X(2)  
      write (*,170) (X(J),J=1,N)    
      write (*,180) (I,F(I),I=1,ME) 
      write (*,190) RES 
      STOP  
      END   


