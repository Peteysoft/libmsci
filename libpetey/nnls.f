       double precision FUNCTION DIFF(X,Y)
c
c  Function used in tests that depend on machine precision.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 7, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
      double precision X, Y
      DIFF=X-Y
      RETURN
      END

      SUBROUTINE G1 (A,B,CTERM,STERM,SIG)   
c
C     COMPUTE ORTHOGONAL ROTATION MATRIX..  
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C   
C     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))   
C                        (-S,C)         (-S,C)(B)   (   0          )    
C     COMPUTE SIG = SQRT(A**2+B**2) 
C        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT 
C        SIG MAY BE IN THE SAME LOCATION AS A OR B .
C     ------------------------------------------------------------------
      double precision A, B, CTERM, ONE, SIG, STERM, XR, YR, ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      if (abs(A) .gt. abs(B)) then
         XR=B/A
         YR=sqrt(ONE+XR**2)
         CTERM=sign(ONE/YR,A)
         STERM=CTERM*XR
         SIG=abs(A)*YR     
         RETURN
      endif

      if (B .ne. ZERO) then
         XR=A/B
         YR=sqrt(ONE+XR**2)
         STERM=sign(ONE/YR,B)
         CTERM=STERM*XR
         SIG=abs(B)*YR     
         RETURN
      endif

      SIG=ZERO  
      CTERM=ZERO  
      STERM=ONE   
      RETURN
      END   


      SUBROUTINE G2    (CTERM,STERM,X,Y)
c
C  APPLY THE ROTATION COMPUTED BY G1 TO (X,Y).  
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1972 DEC 15, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
c     ------------------------------------------------------------------
      double precision CTERM, STERM, X, XR, Y
c     ------------------------------------------------------------------
      XR=CTERM*X+STERM*Y    
      Y=-STERM*X+CTERM*Y    
      X=XR  
      RETURN
      END 
  
            double precision FUNCTION   GEN(ANOISE)
c
C  GENERATE NUMBERS FOR CONSTRUCTION OF TEST CASES. 
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1972 DEC 15, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
c     ------------------------------------------------------------------
      integer I, J, MI, MJ
      double precision AI, AJ, ANOISE, ZERO
      parameter(ZERO = 0.0d0)
      SAVE
c     ------------------------------------------------------------------
      IF (ANOISE) 10,30,20  
   10 MI=891
      MJ=457
      I=5   
      J=7   
      AJ= ZERO
      GEN= ZERO
      RETURN
C   
C     THE SEQUENCE OF VALUES OF J  IS BOUNDED BETWEEN 1 AND 996 
C     IF INITIAL J = 1,2,3,4,5,6,7,8, OR 9, THE PERIOD IS 332   
   20 J=J*MJ
      J=J-997*(J/997)   
      AJ=J-498  
C     THE SEQUENCE OF VALUES OF I  IS BOUNDED BETWEEN 1 AND 999 
C     IF INITIAL I = 1,2,3,6,7, OR 9,  THE PERIOD WILL BE 50
C     IF INITIAL I = 4 OR 8   THE PERIOD WILL BE 25 
C     IF INITIAL I = 5        THE PERIOD WILL BE 10 
   30 I=I*MI
      I=I-1000*(I/1000) 
      AI=I-500  
      GEN=AI+AJ*ANOISE  
      RETURN
      END   


C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
C   
C  CONSTRUCTION AND/OR APPLICATION OF A SINGLE   
C  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B   
C   
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
c                     Subroutine Arguments
c
C     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
c            Householder transformation, or Algorithm H2 to apply a
c            previously constructed transformation.
C     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. 
C     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO   
C            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M     
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
c            vector.  IUE is the storage increment between elements.  
c            On exit when MODE = 1, U() and UP contain quantities
c            defining the vector U of the Householder transformation.
c            on entry with MODE = 2, U() and UP should contain
c            quantities previously computed with MODE = 1.  These will
c            not be modified during the entry with MODE = 2.   
C     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
c            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
c            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
c            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
C     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().  
C     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().  
C     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0  
C            NO OPERATIONS WILL BE DONE ON C(). 
C     ------------------------------------------------------------------
      SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
C     ------------------------------------------------------------------
      integer I, I2, I3, I4, ICE, ICV, INCR, IUE, J
      integer L1, LPIVOT, M, MODE, NCV
      double precision B, C(*), CL, CLINV, ONE, SM
c     double precision U(IUE,M)
      double precision U(IUE,*)
      double precision UP
      parameter(ONE = 1.0d0)
C     ------------------------------------------------------------------
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN    
      CL=abs(U(1,LPIVOT))   
      IF (MODE.EQ.2) GO TO 60   
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
          DO 10 J=L1,M  
   10     CL=MAX(abs(U(1,J)),CL)  
      IF (CL) 130,130,20
   20 CLINV=ONE/CL  
      SM=(U(1,LPIVOT)*CLINV)**2   
          DO 30 J=L1,M  
   30     SM=SM+(U(1,J)*CLINV)**2 
      CL=CL*SQRT(SM)   
      IF (U(1,LPIVOT)) 50,50,40     
   40 CL=-CL
   50 UP=U(1,LPIVOT)-CL 
      U(1,LPIVOT)=CL    
      GO TO 70  
C            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C   
   60 IF (CL) 130,130,70
   70 IF (NCV.LE.0) RETURN  
      B= UP*U(1,LPIVOT)
C                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C   
      IF (B) 80,130,130 
   80 B=ONE/B   
      I2=1-ICV+ICE*(LPIVOT-1)   
      INCR=ICE*(L1-LPIVOT)  
          DO 120 J=1,NCV
          I2=I2+ICV     
          I3=I2+INCR    
          I4=I3 
          SM=C(I2)*UP
              DO 90 I=L1,M  
              SM=SM+C(I3)*U(1,I)
   90         I3=I3+ICE 
          IF (SM) 100,120,100   
  100     SM=SM*B   
          C(I2)=C(I2)+SM*UP
              DO 110 I=L1,M 
              C(I4)=C(I4)+SM*U(1,I)
  110         I4=I4+ICE 
  120     CONTINUE  
  130 RETURN
      END   


      SUBROUTINE LDP (G,MDG,M,N,H,X,XNORM,W,INDEX,MODE)     
c
C  Algorithm LDP: LEAST DISTANCE PROGRAMMING
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1974 MAR 1, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
      integer I, IW, IWDUAL, IY, IZ, J, JF, M, MDG, MODE, N, NP1
c     integer INDEX(M)  
c     double precision G(MDG,N), H(M), X(N), W(*)  
      integer INDEX(*)  
      double precision G(MDG,*), H(*), X(*), W(*)
      double precision DIFF, FAC, ONE, RNORM, XNORM, ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      IF (N.LE.0) GO TO 120 
          DO 10 J=1,N   
   10     X(J)=ZERO     
      XNORM=ZERO
      IF (M.LE.0) GO TO 110 
C   
C     THE DECLARED DIMENSION OF W() MUST BE AT LEAST (N+1)*(M+2)+2*M.   
C   
C      FIRST (N+1)*M LOCS OF W()   =  MATRIX E FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR F FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR Z FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR Y FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR WDUAL FOR PROBLEM NNLS.    
C     COPY G**T INTO FIRST N ROWS AND M COLUMNS OF E.   
C     COPY H**T INTO ROW N+1 OF E.  
C   
      IW=0  
          DO 30 J=1,M   
              DO 20 I=1,N   
              IW=IW+1   
   20         W(IW)=G(J,I)  
          IW=IW+1   
   30     W(IW)=H(J)    
      JF=IW+1   
C                                STORE N ZEROS FOLLOWED BY A ONE INTO F.
          DO 40 I=1,N   
          IW=IW+1   
   40     W(IW)=ZERO    
      W(IW+1)=ONE   
C   
      NP1=N+1   
      IZ=IW+2   
      IY=IZ+NP1 
      IWDUAL=IY+M   
C   
      CALL NNLS (W,NP1,NP1,M,W(JF),W(IY),RNORM,W(IWDUAL),W(IZ),INDEX,   
     *           MODE)  
C                      USE THE FOLLOWING RETURN IF UNSUCCESSFUL IN NNLS.
      IF (MODE.NE.1) RETURN 
      IF (RNORM) 130,130,50 
   50 FAC=ONE   
      IW=IY-1   
          DO 60 I=1,M   
          IW=IW+1   
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
   60     FAC=FAC-H(I)*W(IW)
C   
      IF (DIFF(ONE+FAC,ONE)) 130,130,70 
   70 FAC=ONE/FAC   
          DO 90 J=1,N   
          IW=IY-1   
              DO 80 I=1,M   
              IW=IW+1   
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
   80         X(J)=X(J)+G(I,J)*W(IW)
   90     X(J)=X(J)*FAC 
          DO 100 J=1,N  
  100     XNORM=XNORM+X(J)**2   
      XNORM=sqrt(XNORM) 
C                             SUCCESSFUL RETURN.
  110 MODE=1
      RETURN
C                             ERROR RETURN.       N .LE. 0. 
  120 MODE=2
      RETURN
C                             RETURNING WITH CONSTRAINTS NOT COMPATIBLE.
  130 MODE=4
      RETURN
      END   
 

C     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE)
C   
C  Algorithm NNLS: NONNEGATIVE LEAST SQUARES
C   
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 15, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
c
C     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
C     N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM   
C   
C                      A * X = B  SUBJECT TO X .GE. 0   
C     ------------------------------------------------------------------
c                     Subroutine Arguments
c
C     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE   
C                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N    
C                     MATRIX, A.           ON EXIT A() CONTAINS 
C                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN   
C                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY  
C                     THIS SUBROUTINE.  
C     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- 
C             TAINS Q*B.
C     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL   
C             CONTAIN THE SOLUTION VECTOR.  
C     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE  
C             RESIDUAL VECTOR.  
C     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN    
C             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.  
C             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z   
C     ZZ()     AN M-ARRAY OF WORKING SPACE.     
C     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
C                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS    
C                 P AND Z AS FOLLOWS..  
C   
C                 INDEX(1)   THRU INDEX(NSETP) = SET P.     
C                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.     
C                 IZ1 = NSETP + 1 = NPP1
C                 IZ2 = N   
C     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING 
C             MEANINGS. 
C             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
C             2     THE DIMENSIONS OF THE PROBLEM ARE BAD.  
C                   EITHER M .LE. 0 OR N .LE. 0.
C             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. 
C   
C     ------------------------------------------------------------------
      SUBROUTINE NNLS (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) 
C     ------------------------------------------------------------------
      integer I, II, IP, ITER, ITMAX, IZ, IZ1, IZ2, IZMAX, J, JJ, JZ, L
      integer M, MDA, MODE,N, NPP1, NSETP, RTNKEY
c     integer INDEX(N)  
c     double precision A(MDA,N), B(M), W(N), X(N), ZZ(M)   
      integer INDEX(*)  
      double precision A(MDA,*), B(*), W(*), X(*), ZZ(*)   
      double precision ALPHA, ASAVE, CC, DIFF, DUMMY, FACTOR, RNORM
      double precision SM, SS, T, TEMP, TWO, UNORM, UP, WMAX
      double precision ZERO, ZTEST
      parameter(FACTOR = 0.01d0)
      parameter(TWO = 2.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      MODE=1
      IF (M .le. 0 .or. N .le. 0) then
         MODE=2
         RETURN
      endif
      ITER=0
      ITMAX=3*N 
C   
C                    INITIALIZE THE ARRAYS INDEX() AND X(). 
C   
          DO 20 I=1,N   
          X(I)=ZERO     
   20     INDEX(I)=I    
C   
      IZ2=N 
      IZ1=1 
      NSETP=0   
      NPP1=1
C                             ******  MAIN LOOP BEGINS HERE  ******     
   30 CONTINUE  
C                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
C                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.    
C   
      IF (IZ1 .GT.IZ2.OR.NSETP.GE.M) GO TO 350   
C   
C         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
C   
      DO 50 IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         SM=ZERO   
         DO 40 L=NPP1,M
   40        SM=SM+A(L,J)*B(L)     
         W(J)=SM   
   50 continue
C                                   FIND LARGEST POSITIVE W(J). 
   60 continue
      WMAX=ZERO 
      DO 70 IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         IF (W(J) .gt. WMAX) then
            WMAX=W(J)     
            IZMAX=IZ  
         endif
   70 CONTINUE  
C   
C             IF WMAX .LE. 0. GO TO TERMINATION.
C             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
C   
      IF (WMAX .le. ZERO) go to 350
      IZ=IZMAX  
      J=INDEX(IZ)   
C   
C     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.    
C     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID  
C     NEAR LINEAR DEPENDENCE.   
C   
      ASAVE=A(NPP1,J)   
      CALL H12 (1,NPP1,NPP1+1,M,A(1,J),1,UP,DUMMY,1,1,0)    
      UNORM=ZERO
      IF (NSETP .ne. 0) then
          DO 90 L=1,NSETP   
   90       UNORM=UNORM+A(L,J)**2     
      endif
      UNORM=sqrt(UNORM) 
      IF (DIFF(UNORM+ABS(A(NPP1,J))*FACTOR,UNORM) .gt. ZERO) then
C   
C        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
C        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).    
C   
         DO 120 L=1,M  
  120        ZZ(L)=B(L)    
         CALL H12 (2,NPP1,NPP1+1,M,A(1,J),1,UP,ZZ,1,1,1)   
         ZTEST=ZZ(NPP1)/A(NPP1,J)  
C   
C                                     SEE IF ZTEST IS POSITIVE  
C   
         IF (ZTEST .gt. ZERO) go to 140
      endif
C   
C     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.  
C     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
C     COEFFS AGAIN.     
C   
      A(NPP1,J)=ASAVE   
      W(J)=ZERO 
      GO TO 60  
C   
C     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
C     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER  
C     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN   
C     COL J,  SET W(J)=0.   
C   
  140 continue
      DO 150 L=1,M  
  150    B(L)=ZZ(L)    
C   
      INDEX(IZ)=INDEX(IZ1)  
      INDEX(IZ1)=J  
      IZ1=IZ1+1 
      NSETP=NPP1
      NPP1=NPP1+1   
C   
      IF (IZ1 .le. IZ2) then
         DO 160 JZ=IZ1,IZ2 
            JJ=INDEX(JZ)  
            CALL H12 (2,NSETP,NPP1,M,A(1,J),1,UP,A(1,JJ),1,MDA,1)
  160    continue
      endif
C   
      IF (NSETP .ne. M) then
         DO 180 L=NPP1,M   
  180       A(L,J)=ZERO   
      endif
C   
      W(J)=ZERO 
C                                SOLVE THE TRIANGULAR SYSTEM.   
C                                STORE THE SOLUTION TEMPORARILY IN ZZ().
      RTNKEY = 1
      GO TO 400 
  200 CONTINUE  
C   
C                       ******  SECONDARY LOOP BEGINS HERE ******   
C   
C                          ITERATION COUNTER.   
C 
  210 continue  
      ITER=ITER+1   
      IF (ITER .gt. ITMAX) then
         MODE=3
         write (*,'(/a)') ' NNLS quitting on iteration count.'
         GO TO 350 
      endif
C   
C                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.    
C                                  IF NOT COMPUTE ALPHA.    
C   
      ALPHA=TWO 
      DO 240 IP=1,NSETP 
         L=INDEX(IP)   
         IF (ZZ(IP) .le. ZERO) then
            T=-X(L)/(ZZ(IP)-X(L))     
            IF (ALPHA .gt. T) then
               ALPHA=T   
               JJ=IP 
            endif
         endif
  240 CONTINUE  
C   
C          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL   
C          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.   
C   
      IF (ALPHA.EQ.TWO) GO TO 330   
C   
C          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO   
C          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.    
C   
      DO 250 IP=1,NSETP 
         L=INDEX(IP)   
         X(L)=X(L)+ALPHA*(ZZ(IP)-X(L)) 
  250 continue
C   
C        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I  
C        FROM SET P TO SET Z.   
C   
      I=INDEX(JJ)   
  260 continue
      X(I)=ZERO 
C   
      IF (JJ .ne. NSETP) then
         JJ=JJ+1   
         DO 280 J=JJ,NSETP 
            II=INDEX(J)   
            INDEX(J-1)=II 
            CALL G1 (A(J-1,II),A(J,II),CC,SS,A(J-1,II))   
            A(J,II)=ZERO  
            DO 270 L=1,N  
               IF (L.NE.II) then
c
c                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))  
c
                  TEMP = A(J-1,L)
                  A(J-1,L) = CC*TEMP + SS*A(J,L)
                  A(J,L)   =-SS*TEMP + CC*A(J,L)
               endif
  270       CONTINUE  
c
c                 Apply procedure G2 (CC,SS,B(J-1),B(J))   
c
            TEMP = B(J-1)
            B(J-1) = CC*TEMP + SS*B(J)    
            B(J)   =-SS*TEMP + CC*B(J)    
  280    continue
      endif
c
      NPP1=NSETP
      NSETP=NSETP-1     
      IZ1=IZ1-1 
      INDEX(IZ1)=I  
C   
C        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
C        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
C        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY   
C        THAT ARE NONPOSITIVE WILL BE SET TO ZERO   
C        AND MOVED FROM SET P TO SET Z. 
C   
      DO 300 JJ=1,NSETP 
         I=INDEX(JJ)   
         IF (X(I) .le. ZERO) go to 260
  300 CONTINUE  
C   
C         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
C   
      DO 310 I=1,M  
  310     ZZ(I)=B(I)    
      RTNKEY = 2
      GO TO 400 
  320 CONTINUE  
      GO TO 210 
C                      ******  END OF SECONDARY LOOP  ******
C   
  330 continue
      DO 340 IP=1,NSETP 
          I=INDEX(IP)   
  340     X(I)=ZZ(IP)   
C        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.  
      GO TO 30  
C   
C                        ******  END OF MAIN LOOP  ******   
C   
C                        COME TO HERE FOR TERMINATION.  
C                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.    
C 
  350 continue  
      SM=ZERO   
      IF (NPP1 .le. M) then
         DO 360 I=NPP1,M   
  360       SM=SM+B(I)**2 
      else
         DO 380 J=1,N  
  380       W(J)=ZERO     
      endif
      RNORM=sqrt(SM)    
      RETURN
C   
C     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE     
C     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().     
C   
  400 continue
      DO 430 L=1,NSETP  
         IP=NSETP+1-L  
         IF (L .ne. 1) then
            DO 410 II=1,IP
               ZZ(II)=ZZ(II)-A(II,JJ)*ZZ(IP+1)   
  410       continue
         endif
         JJ=INDEX(IP)  
         ZZ(IP)=ZZ(IP)/A(IP,JJ)    
  430 continue
      go to (200, 320), RTNKEY
      END   

C     SUBROUTINE QRBD (IPASS,Q,E,NN,V,MDV,NRV,C,MDC,NCC)    
c
C  QR ALGORITHM FOR SINGULAR VALUES OF A BIDIAGONAL MATRIX.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------   
C     THE BIDIAGONAL MATRIX 
C   
C                       (Q1,E2,0...    )
C                       (   Q2,E3,0... )
C                D=     (       .      )
C                       (         .   0)
C                       (           .EN)
C                       (          0,QN)
C   
C                 IS PRE AND POST MULTIPLIED BY 
C                 ELEMENTARY ROTATION MATRICES  
C                 RI AND PI SO THAT 
C   
C                 RK...R1*D*P1**(T)...PK**(T) = DIAG(S1,...,SN) 
C   
C                 TO WITHIN WORKING ACCURACY.   
C   
C  1. EI AND QI OCCUPY E(I) AND Q(I) AS INPUT.  
C   
C  2. RM...R1*C REPLACES 'C' IN STORAGE AS OUTPUT.  
C   
C  3. V*P1**(T)...PM**(T) REPLACES 'V' IN STORAGE AS OUTPUT.
C   
C  4. SI OCCUPIES Q(I) AS OUTPUT.   
C   
C  5. THE SI'S ARE NONINCREASING AND NONNEGATIVE.   
C   
C     THIS CODE IS BASED ON THE PAPER AND 'ALGOL' CODE..    
C REF..     
C  1. REINSCH,C.H. AND GOLUB,G.H. 'SINGULAR VALUE DECOMPOSITION 
C     AND LEAST SQUARES SOLUTIONS' (NUMER. MATH.), VOL. 14,(1970).  
C   
C     ------------------------------------------------------------------   
      SUBROUTINE QRBD (IPASS,Q,E,NN,V,MDV,NRV,C,MDC,NCC)    
C     ------------------------------------------------------------------   
      integer MDC, MDV, NCC, NN, NRV
c     double precision C(MDC,NCC), E(NN), Q(NN),V(MDV,NN)
      double precision C(MDC,*  ), E(* ), Q(* ),V(MDV,* )
      integer I, II, IPASS, J, K, KK, L, LL, LP1, N, N10, NQRS
      double precision CS, DIFF, DNORM, F, G, H, SMALL
      double precision ONE, SN, T, TEMP, TWO, X, Y, Z, ZERO
      
      logical WNTV ,HAVERS,FAIL     
      parameter(ONE = 1.0d0, TWO = 2.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------   
      N=NN  
      IPASS=1   
      IF (N.LE.0) RETURN
      N10=10*N  
      WNTV=NRV.GT.0     
      HAVERS=NCC.GT.0   
      FAIL=.FALSE.  
      NQRS=0
      E(1)=ZERO 
      DNORM=ZERO
           DO 10 J=1,N  
   10      DNORM=max(abs(Q(J))+abs(E(J)),DNORM)   
           DO 200 KK=1,N
           K=N+1-KK     
C   
C     TEST FOR SPLITTING OR RANK DEFICIENCIES.. 
C         FIRST MAKE TEST FOR LAST DIAGONAL TERM, Q(K), BEING SMALL.    
   20       IF(K.EQ.1) GO TO 50     
            IF(DIFF(DNORM+Q(K),DNORM) .ne. ZERO) go to 50
C   
C     SINCE Q(K) IS SMALL WE WILL MAKE A SPECIAL PASS TO    
C     TRANSFORM E(K) TO ZERO.   
C   
           CS=ZERO  
           SN=-ONE  
                DO 40 II=2,K
                I=K+1-II
                F=-SN*E(I+1)
                E(I+1)=CS*E(I+1)    
                CALL G1 (Q(I),F,CS,SN,Q(I))     
C         TRANSFORMATION CONSTRUCTED TO ZERO POSITION (I,K).
C   
                IF (.NOT.WNTV) GO TO 40 
                     DO 30 J=1,NRV  
c
c                          Apply procedure G2 (CS,SN,V(J,I),V(J,K))  
c
                        TEMP = V(J,I)
                        V(J,I) = CS*TEMP + SN*V(J,K)
                        V(J,K) =-SN*TEMP + CS*V(J,K)
   30                continue
C              ACCUMULATE RT. TRANSFORMATIONS IN V. 
C   
   40           CONTINUE
C   
C         THE MATRIX IS NOW BIDIAGONAL, AND OF LOWER ORDER  
C         SINCE E(K) .EQ. ZERO..    
C   
   50           DO 60 LL=1,K
                  L=K+1-LL
                  IF(DIFF(DNORM+E(L),DNORM) .eq. ZERO) go to 100
                  IF(DIFF(DNORM+Q(L-1),DNORM) .eq. ZERO) go to 70
   60           CONTINUE
C     THIS LOOP CAN'T COMPLETE SINCE E(1) = ZERO.   
C   
           GO TO 100    
C   
C         CANCELLATION OF E(L), L.GT.1. 
   70      CS=ZERO  
           SN=-ONE  
                DO 90 I=L,K 
                F=-SN*E(I)  
                E(I)=CS*E(I)
                IF(DIFF(DNORM+F,DNORM) .eq. ZERO) go to 100
                CALL G1 (Q(I),F,CS,SN,Q(I))     
                IF (HAVERS) then
                     DO 80 J=1,NCC  
c
c                          Apply procedure G2 ( CS, SN, C(I,J), C(L-1,J)
c
                        TEMP = C(I,J)
                        C(I,J)   = CS*TEMP + SN*C(L-1,J)
                        C(L-1,J) =-SN*TEMP + CS*C(L-1,J)
   80                continue
                endif
   90           CONTINUE
C   
C         TEST FOR CONVERGENCE..    
  100      Z=Q(K)   
           IF (L.EQ.K) GO TO 170    
C   
C         SHIFT FROM BOTTOM 2 BY 2 MINOR OF B**(T)*B.   
           X=Q(L)   
           Y=Q(K-1)     
           G=E(K-1)     
           H=E(K)   
           F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
           G=sqrt(ONE+F**2) 
           IF (F .ge. ZERO) then
              T=F+G
           else
              T=F-G
           endif
           F=((X-Z)*(X+Z)+H*(Y/T-H))/X  
C   
C         NEXT QR SWEEP..   
           CS=ONE   
           SN=ONE   
           LP1=L+1  
                DO 160 I=LP1,K  
                G=E(I)  
                Y=Q(I)  
                H=SN*G  
                G=CS*G  
                CALL G1 (F,H,CS,SN,E(I-1))  
                F=X*CS+G*SN 
                G=-X*SN+G*CS
                H=Y*SN  
                Y=Y*CS  
                IF (WNTV) then
C   
C              ACCUMULATE ROTATIONS (FROM THE RIGHT) IN 'V' 
c
                     DO 130 J=1,NRV 
c                    
c                          Apply procedure G2 (CS,SN,V(J,I-1),V(J,I))
c
                        TEMP = V(J,I-1)
                        V(J,I-1) = CS*TEMP + SN*V(J,I)
                        V(J,I)   =-SN*TEMP + CS*V(J,I)
  130                continue
                endif
                CALL G1 (F,H,CS,SN,Q(I-1))  
                F=CS*G+SN*Y 
                X=-SN*G+CS*Y
                IF (HAVERS) then
                     DO 150 J=1,NCC 
c
c                          Apply procedure G2 (CS,SN,C(I-1,J),C(I,J))
c
                        TEMP = C(I-1,J)
                        C(I-1,J) = CS*TEMP + SN*C(I,J)
                        C(I,J)   =-SN*TEMP + CS*C(I,J)
  150                continue
                endif
c
C              APPLY ROTATIONS FROM THE LEFT TO 
C              RIGHT HAND SIDES IN 'C'..
C   
  160           CONTINUE
           E(L)=ZERO    
           E(K)=F   
           Q(K)=X   
           NQRS=NQRS+1  
           IF (NQRS.LE.N10) GO TO 20
C          RETURN TO 'TEST FOR SPLITTING'.  
C 
           SMALL=ABS(E(K))
           I=K     
C          IF FAILURE TO CONVERGE SET SMALLEST MAGNITUDE
C          TERM IN OFF-DIAGONAL TO ZERO.  CONTINUE ON.
C      ..           
                DO 165 J=L,K 
                TEMP=ABS(E(J))
                IF(TEMP .EQ. ZERO) GO TO 165
                IF(TEMP .LT. SMALL) THEN
                     SMALL=TEMP
                     I=J
                end if  
  165           CONTINUE
           E(I)=ZERO
           NQRS=0
           FAIL=.TRUE.  
           GO TO 20
C     ..    
C     CUTOFF FOR CONVERGENCE FAILURE. 'NQRS' WILL BE 2*N USUALLY.   
  170      IF (Z.GE.ZERO) GO TO 190 
           Q(K)=-Z  
           IF (WNTV) then
                DO 180 J=1,NRV  
  180           V(J,K)=-V(J,K)  
           endif
  190      CONTINUE     
C         CONVERGENCE. Q(K) IS MADE NONNEGATIVE..   
C   
  200      CONTINUE     
      IF (N.EQ.1) RETURN
           DO 210 I=2,N 
           IF (Q(I).GT.Q(I-1)) GO TO 220
  210      CONTINUE     
      IF (FAIL) IPASS=2 
      RETURN
C     ..    
C     EVERY SINGULAR VALUE IS IN ORDER..
  220      DO 270 I=2,N 
           T=Q(I-1)     
           K=I-1
                DO 230 J=I,N
                IF (T.GE.Q(J)) GO TO 230
                T=Q(J)  
                K=J     
  230           CONTINUE
           IF (K.EQ.I-1) GO TO 270  
           Q(K)=Q(I-1)  
           Q(I-1)=T     
           IF (HAVERS) then
                DO 240 J=1,NCC  
                T=C(I-1,J)  
                C(I-1,J)=C(K,J)     
  240           C(K,J)=T
           endif

  250      IF (WNTV) then
                DO 260 J=1,NRV  
                T=V(J,I-1)  
                V(J,I-1)=V(J,K)     
  260           V(J,K)=T
           endif
  270      CONTINUE     
C         END OF ORDERING ALGORITHM.
C   
      IF (FAIL) IPASS=2 
      RETURN
      END   

C     SUBROUTINE SVDRS (A, MDA, M1, N1, B, MDB, NB, S, WORK) 
c
C  SINGULAR VALUE DECOMPOSITION ALSO TREATING RIGHT SIDE VECTOR.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1974 SEP 25, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------   
c  This 1995 version differs from the original 1974 version by adding
c  the argument WORK().
c  WORK() provides 2*N1 locations of work space.  Originally S() was
c  required to have 3*N1 elements, of which the last 2*N1 were used for
c  work space.  Now S() only needs N1 elements.
C     ------------------------------------------------------------------
c  This subroutine computes the singular value decomposition of the
c  given M1 x N1 matrix, A, and optionally applys the transformations
c  from the left to the NB column vectors of the M1 x NB matrix B.
c  Either M1 .ge. N1  or  M1 .lt. N1 is permitted.
c
c       The singular value decomposition of A is of the form
c
c                  A  =  U * S * V**t
c
c  where U is M1 x M1 orthogonal, S is M1 x N1 diagonal with the
c  diagonal terms nonnegative and ordered from large to small, and
c  V is N1 x N1 orthogonal.  Note that these matrices also satisfy
c
c                  S = (U**t) * A * V
c
c       The matrix V is returned in the leading N1 rows and
c  columns of the array A(,).
c
c       The singular values, i.e. the diagonal terms of the matrix S,
c  are returned in the array S().  If M1 .lt. N1, positions M1+1
c  through N1 of S() will be set to zero.
c
c       The product matrix  G = U**t * B replaces the given matrix B
c  in the array B(,).
c
c       If the user wishes to obtain a minimum length least squares
c  solution of the linear system
c
c                        A * X ~=~ B
c
c  the solution X can be constructed, following use of this subroutine,
c  by computing the sum for i = 1, ..., R of the outer products
c
c          (Col i of V) * (1/S(i)) * (Row i of G)
c
c  Here R denotes the pseudorank of A which the user may choose
c  in the range 0 through Min(M1, N1) based on the sizes of the
c  singular values.
C     ------------------------------------------------------------------
C                    Subroutine Arguments
c
c  A(,)     (In/Out)  On input contains the M1 x N1 matrix A.
c           On output contains the N1 x N1 matrix V.
c
c  LDA      (In)  First dimensioning parameter for A(,).
c           Require LDA .ge. Max(M1, N1).
c
c  M1       (In)  No. of rows of matrices A, B, and G.
c           Require M1 > 0.
c
c  N1       (In)  No. of cols of matrix A, No. of rows and cols of
c           matrix V.  Permit M1 .ge. N1  or  M1 .lt. N1.
c           Require N1 > 0.
c
c  B(,)     (In/Out)  If NB .gt. 0 this array must contain an
c           M1 x NB matrix on input and will contain the
c           M1 x NB product matrix, G = (U**t) * B on output.
c
c  LDB      (In)  First dimensioning parameter for B(,).
c           Require LDB .ge. M1.
c
c  NB       (In)  No. of cols in the matrices B and G.
c           Require NB .ge. 0.
c
c  S()      (Out)  Must be dimensioned at least N1.  On return will
c           contain the singular values of A, with the ordering
c                S(1) .ge. S(2) .ge. ... .ge. S(N1) .ge. 0.
c           If M1 .lt. N1 the singular values indexed from M1+1
c           through N1 will be zero.
c           If the given integer arguments are not consistent, this
c           subroutine will return immediately, setting S(1) = -1.0.
c
c  WORK()  (Scratch)  Work space of total size at least 2*N1.
c           Locations 1 thru N1 will hold the off-diagonal terms of
c           the bidiagonal matrix for subroutine QRBD.  Locations N1+1
c           thru 2*N1 will save info from one call to the next of
c           H12.
c     ------------------------------------------------------------------
C  This code gives special treatment to rows and columns that are
c  entirely zero.  This causes certain zero sing. vals. to appear as
c  exact zeros rather than as about MACHEPS times the largest sing. val.
c  It similarly cleans up the associated columns of U and V.  
c
c  METHOD..  
c   1. EXCHANGE COLS OF A TO PACK NONZERO COLS TO THE LEFT.  
c      SET N = NO. OF NONZERO COLS.  
c      USE LOCATIONS A(1,N1),A(1,N1-1),...,A(1,N+1) TO RECORD THE    
c      COL PERMUTATIONS. 
c   2. EXCHANGE ROWS OF A TO PACK NONZERO ROWS TO THE TOP.   
c      QUIT PACKING IF FIND N NONZERO ROWS.  MAKE SAME ROW EXCHANGES 
c      IN B.  SET M SO THAT ALL NONZERO ROWS OF THE PERMUTED A   
c      ARE IN FIRST M ROWS.  IF M .LE. N THEN ALL M ROWS ARE 
c      NONZERO.  IF M .GT. N THEN      THE FIRST N ROWS ARE KNOWN    
c      TO BE NONZERO,AND ROWS N+1 THRU M MAY BE ZERO OR NONZERO.     
c   3. APPLY ORIGINAL ALGORITHM TO THE M BY N PROBLEM.   
c   4. MOVE PERMUTATION RECORD FROM A(,) TO S(I),I=N+1,...,N1.   
c   5. BUILD V UP FROM  N BY N  TO  N1 BY N1  BY PLACING ONES ON     
c      THE DIAGONAL AND ZEROS ELSEWHERE.  THIS IS ONLY PARTLY DONE   
c      EXPLICITLY.  IT IS COMPLETED DURING STEP 6.   
c   6. EXCHANGE ROWS OF V TO COMPENSATE FOR COL EXCHANGES OF STEP 2. 
c   7. PLACE ZEROS IN  S(I),I=N+1,N1  TO REPRESENT ZERO SING VALS.   
c     ------------------------------------------------------------------
      subroutine SVDRS (A, MDA, M1, N1, B, MDB, NB, S, WORK) 
      integer I, IPASS, J, K, L, M, MDA, MDB, M1
      integer N, NB, N1, NP1, NS, NSP1
c     double precision A(MDA,N1),B(MDB,NB), S(N1)
      double precision A(MDA, *),B(MDB, *), S( *)
      double precision ONE, T, WORK(N1,2), ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                             BEGIN.. SPECIAL FOR ZERO ROWS AND COLS.   
C   
C                             PACK THE NONZERO COLS TO THE LEFT 
C   
      N=N1  
      IF (N.LE.0.OR.M1.LE.0) RETURN 
      J=N   
   10 CONTINUE  
         DO 20 I=1,M1  
            IF (A(I,J) .ne. ZERO) go to 50
   20    CONTINUE  
C   
C           COL J  IS ZERO. EXCHANGE IT WITH COL N.   
C   
         IF (J .ne. N) then
            DO 30 I=1,M1  
   30       A(I,J)=A(I,N) 
         endif
         A(1,N)=J  
         N=N-1 
   50    CONTINUE  
         J=J-1 
      IF (J.GE.1) GO TO 10  
C                             IF N=0 THEN A IS ENTIRELY ZERO AND SVD    
C                             COMPUTATION CAN BE SKIPPED    
      NS=0  
      IF (N.EQ.0) GO TO 240 
C                             PACK NONZERO ROWS TO THE TOP  
C                             QUIT PACKING IF FIND N NONZERO ROWS   
      I=1   
      M=M1  
   60 IF (I.GT.N.OR.I.GE.M) GO TO 150   
      IF (A(I,I)) 90,70,90  
   70     DO 80 J=1,N   
          IF (A(I,J)) 90,80,90  
   80     CONTINUE  
      GO TO 100 
   90 I=I+1 
      GO TO 60  
C                             ROW I IS ZERO     
C                             EXCHANGE ROWS I AND M 
  100 IF(NB.LE.0) GO TO 115 
          DO 110 J=1,NB 
          T=B(I,J)  
          B(I,J)=B(M,J) 
  110     B(M,J)=T  
  115     DO 120 J=1,N  
  120     A(I,J)=A(M,J) 
      IF (M.GT.N) GO TO 140 
          DO 130 J=1,N  
  130     A(M,J)=ZERO   
  140 CONTINUE  
C                             EXCHANGE IS FINISHED  
      M=M-1 
      GO TO 60  
C   
  150 CONTINUE  
C                             END.. SPECIAL FOR ZERO ROWS AND COLUMNS   
C                             BEGIN.. SVD ALGORITHM 
C     METHOD..  
C     (1)     REDUCE THE MATRIX TO UPPER BIDIAGONAL FORM WITH   
C     HOUSEHOLDER TRANSFORMATIONS.  
C          H(N)...H(1)AQ(1)...Q(N-2) = (D**T,0)**T  
C     WHERE D IS UPPER BIDIAGONAL.  
C   
C     (2)     APPLY H(N)...H(1) TO B.  HERE H(N)...H(1)*B REPLACES B    
C     IN STORAGE.   
C   
C     (3)     THE MATRIX PRODUCT W= Q(1)...Q(N-2) OVERWRITES THE FIRST  
C     N ROWS OF A IN STORAGE.   
C   
C     (4)     AN SVD FOR D IS COMPUTED.  HERE K ROTATIONS RI AND PI ARE 
C     COMPUTED SO THAT  
C          RK...R1*D*P1**(T)...PK**(T) = DIAG(S1,...,SM)    
C     TO WORKING ACCURACY.  THE SI ARE NONNEGATIVE AND NONINCREASING.   
C     HERE RK...R1*B OVERWRITES B IN STORAGE WHILE  
C     A*P1**(T)...PK**(T)  OVERWRITES A IN STORAGE. 
C   
C     (5)     IT FOLLOWS THAT,WITH THE PROPER DEFINITIONS,  
C     U**(T)*B OVERWRITES B, WHILE V OVERWRITES THE FIRST N ROW AND     
C     COLUMNS OF A.     
C   
      L=min(M,N)   
C             THE FOLLOWING LOOP REDUCES A TO UPPER BIDIAGONAL AND  
C             ALSO APPLIES THE PREMULTIPLYING TRANSFORMATIONS TO B.     
C   
          DO 170 J=1,L  
          IF (J.GE.M) GO TO 160     
          CALL H12 (1,J,J+1,M,A(1,J),1,T,A(1,J+1),1,MDA,N-J)
          CALL H12 (2,J,J+1,M,A(1,J),1,T,B,1,MDB,NB)
  160     IF (J.GE.N-1) GO TO 170   
          CALL H12 (1,J+1,J+2,N,A(J,1),MDA,work(J,2),A(J+1,1),MDA,1,M-J)   
  170     CONTINUE  
C   
C     COPY THE BIDIAGONAL MATRIX INTO S() and WORK() FOR QRBD.   
C 1986 Jan 8. C. L. Lawson. Changed N to L in following 2 statements.
      IF (L.EQ.1) GO TO 190 
          DO 180 J=2,L  
          S(J)=A(J,J) 
  180     WORK(J,1)=A(J-1,J)   
  190 S(1)=A(1,1)     
C   
      NS=N  
      IF (M.GE.N) GO TO 200 
      NS=M+1
      S(NS)=ZERO  
      WORK(NS,1)=A(M,M+1)  
  200 CONTINUE  
C   
C     CONSTRUCT THE EXPLICIT N BY N PRODUCT MATRIX, W=Q1*Q2*...*QL*I    
C     IN THE ARRAY A(). 
C   
          DO 230 K=1,N  
          I=N+1-K   
          IF (I .GT. min(M,N-2)) GO TO 210     
          CALL H12 (2,I+1,I+2,N,A(I,1),MDA,WORK(I,2),A(1,I+1),1,MDA,N-I)   
  210         DO 220 J=1,N  
  220         A(I,J)=ZERO   
  230     A(I,I)=ONE    
C   
C          COMPUTE THE SVD OF THE BIDIAGONAL MATRIX 
C   
      CALL QRBD (IPASS,S(1),WORK(1,1),NS,A,MDA,N,B,MDB,NB)   
C   
      if(IPASS .eq. 2) then
         write (*,'(/a)')
     *      ' FULL ACCURACY NOT ATTAINED IN BIDIAGONAL SVD'
      endif

  240 CONTINUE  
      IF (NS.GE.N) GO TO 260
      NSP1=NS+1 
          DO 250 J=NSP1,N   
  250     S(J)=ZERO   
  260 CONTINUE  
      IF (N.EQ.N1) RETURN   
      NP1=N+1   
C                                  MOVE RECORD OF PERMUTATIONS  
C                                  AND STORE ZEROS  
          DO 280 J=NP1,N1   
          S(J)=A(1,J) 
              DO 270 I=1,N  
  270         A(I,J)=ZERO   
  280     CONTINUE  
C                             PERMUTE ROWS AND SET ZERO SINGULAR VALUES.
          DO 300 K=NP1,N1   
          I=S(K)  
          S(K)=ZERO   
              DO 290 J=1,N1 
              A(K,J)=A(I,J) 
  290         A(I,J)=ZERO   
          A(I,K)=ONE    
  300     CONTINUE  
C                             END.. SPECIAL FOR ZERO ROWS AND COLUMNS   
      RETURN
      END   


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


