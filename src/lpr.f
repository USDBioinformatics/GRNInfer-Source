C     ALGORITHM 478 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN COMM. ACM, VOL. 17, NO. 06,
C     P. 319.
      SUBROUTINE L1(M,N,M2,N2,A,B,TOLER,X,E,S)                             A 010
C THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX METHOD
C OF LINEAR PROGRAMMING TO CALCULATE AN L1 SOLUTION TO AN
C OVER-DETERMINED SYSTEM OF LINEAR EQUATIONS.
C DESCRIPTION OF PARAMETERS.
C M      NUMBER OF EQUATIONS.
C N      NUMBER OF UNKNOWNS (M.GE.N).
C M2     SET EQUAL TO M+2 FOR ADJUSTABLE DIMENSIONS.
C N2     SET EQUAL TO N+2 FOR ADJUSTABLE DIMENSIONS.
C A      TWO DIMENSIONAL REAL ARRAY OF SIZE (M2,N2).
C        ON ENTRY, THE COEFFICIENTS OF THE MATRIX MUST BE
C        STORED IN THE FIRST M ROWS AND N COLUMNS OF A.
C        THESE VALUES ARE DESTROYED BY THE SUBROUTINE.
C B      ONE DIMENSIONAL REAL ARRAY OF SIZE M. ON ENTRY, B
C        MUST CONTAIN THE RIGHT HAND SIDE OF THE EQUATIONS.
C        THESE VALUES ARE DESTROYED BY THE SUBROUTINE.
C TOLER  A SMALL POSITIVE TOLERANCE. EMPIRICAL EVIDENCE
C        SUGGESTS TOLER=10**(-D*2/3) WHERE D REPRESENTS
C        THE NUMBER OF DECIMAL DIGITS OF ACCURACY AVALABLE
C        (SEE DESCRIPTION).
C X      ONE DIMENSIONAL REAL ARRAY OF SIZE N. ON EXIT, THIS
C        ARRAY CONTAINS A SOLUTION TO THE L1 PROBLEM.
C E      ONE DIMENSIONAL REAL ARRAY OF SIZE M. ON EXIT, THIS
C        ARRAY CONTAINS THE RESIDUALS IN THE EQUATIONS.
C S      INTEGER ARRAY OF SIZE M USED FOR WORKSPACE.
C ON EXIT FROM THE SUBROUTINE, THE ARRAY A CONTAINS THE
C FOLLOWING INFORMATION.
C A(M+1,N+1)  THE MINIMUM SUM OF THE ABSOLUTE VALUES OF
C             THE RESIDUALS.
C A(M+1,N+2)  THE RANK OF THE MATRIX OF COEFFICIENTS.
C A(M+2,N+1)  EXIT CODE WITH VALUES.
C             0 - OPTIMAL SOLUTION WHICH IS PROBABLY NON-
C                 UNIQUE (SEE DESCRIPTION).
C             1 - UNIQUE OPTIMAL SOLUTION.
C             2 - CALCULATIONS TERMINATED PREMATURELY DUE TO
C                 ROUNDING ERRORS.
C A(M+2,N+2)  NUMBER OF SIMPLEX ITERATIONS PERFORMED.
      DOUBLE PRECISION SUM
      DOUBLE PRECISION MIN, MAX, A(M2,N2), X(N), E(M), B(M)
      INTEGER OUT, S(M)
      LOGICAL STAGE, TEST
C BIG MUST BE SET EQUAL TO ANY VERY LARGE REAL CONSTANT.
C ITS VALUE HERE IS APPROPRIATE FOR THE IBM 370.
      DATA BIG/1.E38/
C INITIALIZATION.
      M1 = M + 1
      N1 = N + 1
      DO 10 J=1,N
        A(M2,J) = J
        X(J) = 0.
   10 CONTINUE
      DO 40 I=1,M
        A(I,N2) = N + I
        A(I,N1) = B(I)
        IF (B(I).GE.0.) GO TO 30
        DO 20 J=1,N2
          A(I,J) = -A(I,J)
   20   CONTINUE
   30   E(I) = 0.
   40 CONTINUE
C COMPUTE THE MARGINAL COSTS.
      DO 60 J=1,N1
        SUM = 0.D0
        DO 50 I=1,M
          SUM = SUM + A(I,J)
   50   CONTINUE
        A(M1,J) = SUM
   60 CONTINUE
C STAGE I.
C	print *, 'Stage I'
C DETERMINE THE VECTOR TO ENTER THE BASIS.
      STAGE = .TRUE.
      KOUNT = 0
      KR = 1
      KL = 1
   70 MAX = -1.
      DO 80 J=KR,N
        IF (ABS(A(M2,J)).GT.N) GO TO 80
        D = ABS(A(M1,J))
        IF (D.LE.MAX) GO TO 80
        MAX = D
        IN = J
   80 CONTINUE
      IF (A(M1,IN).GE.0.) GO TO 100
      DO 90 I=1,M2
        A(I,IN) = -A(I,IN)
   90 CONTINUE
C DETERMINE THE VECTOR TO LEAVE THE BASIS.
  100 K = 0
      DO 110 I=KL,M
        D = A(I,IN)
        IF (D.LE.TOLER) GO TO 110
        K = K + 1
        B(K) = A(I,N1)/D
        S(K) = I
        TEST = .TRUE.
  110 CONTINUE
  120 IF (K.GT.0) GO TO 130
      TEST = .FALSE.
      GO TO 150
  130 MIN = BIG
      DO 140 I=1,K
        IF (B(I).GE.MIN) GO TO 140
        J = I
        MIN = B(I)
        OUT = S(I)
  140 CONTINUE
      B(J) = B(K)
      S(J) = S(K)
      K = K - 1
C CHECK FOR LINEAR DEPENDENCE IN STAGE I.
  150 IF (TEST .OR. .NOT.STAGE) GO TO 170
      DO 160 I=1,M2
        D = A(I,KR)
        A(I,KR) = A(I,IN)
        A(I,IN) = D
  160 CONTINUE
      KR = KR + 1
      GO TO 260
  170 IF (TEST) GO TO 180
      A(M2,N1) = 2.
      GO TO 350
  180 PIVOT = A(OUT,IN)
      IF (A(M1,IN)-PIVOT-PIVOT.LE.TOLER) GO TO 200
      DO 190 J=KR,N1
        D = A(OUT,J)
        A(M1,J) = A(M1,J) - D - D
        A(OUT,J) = -D
  190 CONTINUE
      A(OUT,N2) = -A(OUT,N2)
      GO TO 120
C PIVOT ON A(OUT,IN).
  200 DO 210 J=KR,N1
        IF (J.EQ.IN) GO TO 210
        A(OUT,J) = A(OUT,J)/PIVOT
  210 CONTINUE
      DO 230 I=1,M1
        IF (I.EQ.OUT) GO TO 230
        D = A(I,IN)
        DO 220 J=KR,N1
          IF (J.EQ.IN) GO TO 220
          A(I,J) = A(I,J) - D*A(OUT,J)
  220   CONTINUE
  230 CONTINUE
      DO 240 I=1,M1
        IF (I.EQ.OUT) GO TO 240
        A(I,IN) = -A(I,IN)/PIVOT
  240 CONTINUE
      A(OUT,IN) = 1./PIVOT
      D = A(OUT,N2)
      A(OUT,N2) = A(M2,IN)
      A(M2,IN) = D
      KOUNT = KOUNT + 1
      IF (.NOT.STAGE) GO TO 270
C INTERCHANGE ROWS IN STAGE I.
      KL = KL + 1
      DO 250 J=KR,N2
        D = A(OUT,J)
        A(OUT,J) = A(KOUNT,J)
        A(KOUNT,J) = D
  250 CONTINUE
  260 IF (KOUNT+KR.NE.N1) GO TO 70
C STAGE II.
C	print *, 'Stage II'
      STAGE = .FALSE.
C DETERMINE THE VECTOR TO ENTER THE BASIS.
  270 MAX = -BIG
      DO 290 J=KR,N
        D = A(M1,J)
        IF (D.GE.0.) GO TO 280
        IF (D.GT.(-2.)) GO TO 290
        D = -D - 2.
  280   IF (D.LE.MAX) GO TO 290
        MAX = D
        IN = J
  290 CONTINUE
      IF (MAX.LE.TOLER) GO TO 310
      IF (A(M1,IN).GT.0.) GO TO 100
      DO 300 I=1,M2
        A(I,IN) = -A(I,IN)
  300 CONTINUE
      A(M1,IN) = A(M1,IN) - 2.
      GO TO 100
C PREPARE OUTPUT.
	print *, 'Prepare Output'
  310 L = KL - 1
      DO 330 I=1,L
        IF (A(I,N1).GE.0.) GO TO 330
        DO 320 J=KR,N2
          A(I,J) = -A(I,J)
  320   CONTINUE
  330 CONTINUE
      A(M2,N1) = 0.
      IF (KR.NE.1) GO TO 350
      DO 340 J=1,N
        D = ABS(A(M1,J))
        IF (D.LE.TOLER .OR. 2.-D.LE.TOLER) GO TO 350
  340 CONTINUE
      A(M2,N1) = 1.
  350 DO 380 I=1,M
        K = A(I,N2)
        D = A(I,N1)
        IF (K.GT.0) GO TO 360
        K = -K
        D = -D
  360   IF (I.GE.KL) GO TO 370
        X(K) = D
        GO TO 380
  370   K = K - N
        E(K) = D
  380 CONTINUE
      A(M2,N2) = KOUNT
      A(M1,N2) = N1 - KR
      SUM = 0.D0
      DO 390 I=KL,M
        SUM = SUM + A(I,N1)
  390 CONTINUE
      A(M1,N1) = SUM
      RETURN
      END
