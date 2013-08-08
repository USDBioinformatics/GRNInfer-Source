!  regnw.f90 
!
!  FUNCTIONS:
!	regnw      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: regnw
!
!  PURPOSE:  Reverse Engineering gene networks by combining multiple datasets.
!  Author    Yong Wang    2005/12
!****************************************************************************
PROGRAM regnw 
	CHARACTER(100) :: command_arg
	Integer :: command_arg_ct
	Integer :: i

	command_arg_ct = COMMAND_ARGUMENT_COUNT()
	
	DO i=1,command_arg_ct
		CALL GET_COMMAND_ARGUMENT(i,command_arg)
		call mainprocess(TRIM(command_arg))
	END DO
END PROGRAM regnw

SUBROUTINE mainprocess(inputFileName)
	USE nrtype
	USE nr
	USE large_array

	IMPLICIT NONE

	!	The first part: reading the data and performing the SVD decomposition
	INTEGER(I4B), PARAMETER :: MP=1000, NP=66, ND=10   !Max of Number of protein, samples and dataset respectively
	INTEGER(I4B) :: i,j,m,n,k,q,dc,NumDs, Location(1), Nzero,IteNum,IteMax, RealorSim
	real(DP) :: NumSapAll, lamda, epsilon,Toler, ObjFun, PreObjFun,temp, SubObj, SubMin,SumOmega, AveV, AveSigma, Threshold !,SubObjTemp,SubMinTemp,Coeff
	INTEGER(I4B), DIMENSION(ND) :: NumSample, L
	Double Precision, DIMENSION(ND) :: Weight, Wdata
	CHARACTER(3) :: dummy
	CHARACTER(LEN=*), INTENT(IN) :: inputFileName
	CHARACTER(LEN=LEN(inputFileName)-4) :: inputFile

	inputFile=inputFileName

	lamda=0.0               !The parameter to achieve the balance between sparse and objective function
	epsilon=0.000001d0      !The parameter to control the precision of computation
	Threshold= 0.000001d0   !The parameter is precision of the inferred edge strength
	IteMax=100              !The paramter to control the maximal iteration number
			           
	RealorSim=1             ! 1 denote the real dataset 2**x, 0 denote simulated data

	OPEN(UNIT=7, FILE=inputFile//'.txt',STATUS='old')

	OPEN(UNIT=8, FILE=inputFile//'output.dat',STATUS='REPLACE')

	READ(7,*)
	READ(7,*) NumDs
	READ(7,*)
	READ(7,*) m
	READ(7,*)

	!Here is the right place to allocate the JK and v
	ALLOCATE(GeneName(m))
	ALLOCATE(JK(NumDs,m,m))
	ALLOCATE(v(NumDs,m,m))
	ALLOCATE(GNW(m,m))
	JK=0.0
	v=0.0
	GNW=0.0

	DO i=1,m
		READ(7,*)  GeneName(i)
	END DO

	NumDs=0

	DO
		READ(7,'(a)') dummy
		IF (dummy == 'END') EXIT
		NumDs=NumDs+1
		READ(7,*)          !Fro the brief introduction of dataset
		READ(7,*)          !For the dimension of the microarray matrix
		READ(7,*) m, NumSample(NumDs)
		n=NumSample(NumDs)-1

		!When the m and n is prepared it is time for allocating the array used by svd
		ALLOCATE(a(m,n+1))
		ALLOCATE(B(m,n))
		ALLOCATE(DifX(m,n))  
		ALLOCATE(Ts(n+1))
		ALLOCATE(w(m))
		ALLOCATE(Uw(m))
		ALLOCATE(u(n+1,m))
		ALLOCATE(A0(n,m))
		a=0.0
		B=0.0
		DifX=0.0
		Ts=0.0
		w=0.0
		Uw=0.0
		u=0.0
		A0=0.0		
		
		READ(7,*)              !For the time point of every sample
		READ(7,*) (Ts(j), j=1,NumSample(NumDs))
		READ(7,*)              !For the weight of every dataset
		READ(7,*) Wdata(NumDs)
		READ(7,*)              !For the matrix seciton 

		!copy original matrix into u
		DO i=1,m
			READ(7,*) (a(i,j), j=1,NumSample(NumDs))
			IF (RealorSim .EQ. 1) THEN
				u(1:NumSample(NumDs), i)=2**a(i,1:NumSample(NumDs))    !This is for real data 2**a(NumDs,i,1:NumSample(NumDs))  
			ELSE
				u(1:NumSample(NumDs), i)=a(i,1:NumSample(NumDs))       !This is for simulated data
			END IF
		END DO

		B(1:m,1:n)=0.0
		DO j=1, n
			DifX(1:m,j)=(u(j+1,1:m)-u(j,1:m))/(Ts(j+1)-Ts(j))           !Record the X matrix before changed in svdcmp
		END DO

		DO i=1,n
			WRITE(8,'(1x,6f12.6)') (u(i, j),j=1,m)
		END DO

		!perform decomposition
		CALL dsvdcmp(u,n,m,n+1,m,w,GNW)           !Here GNW as a temp variable to solve the stack flow problem
		v(NumDs,1:m,1:m)=GNW

		!Sorting of singular value need here
		i=1
		temp=maxval(w(1:m))
		DO WHILE ((temp>0.0) .AND. (i<=n))
			Location=MAXLOC(w(1:m-i+1)) 

			!swap w
			temp=w(m-i+1)
			w(m-i+1) = w(Location(1))
			w(Location(1))=temp

			!swap u
			DO j=1,n
				temp=u(j,m-i+1)
				u(j, m-i+1) = u(j,Location(1))
				u(j, Location(1))=temp
			END DO

			!swap v
			DO 	j=1,m
				temp=v(NumDs, j, m-i+1)
				v(NumDs, j,m-i+1) = v(NumDs, j, Location(1))
				v(NumDs,j, Location(1))=temp
			END DO
			i=i+1
			temp=maxval(w(1:m-i))
		END DO


		!construct the particular solution of every dataset by svd results
		L(NumDs)= 0   !m-n
		DO i=1,m
			IF (w(i)<=0.000001) THEN      !set the threshold of svd is epsilon
				Uw(i)=0.0
				L(NumDs)=L(NumDs)+1
			ELSE 
			    Uw(i)=1.0/w(i)
			END IF
		END DO
		
		NumSapAll=NumSapAll+NumSample(NumDs)

		DEALLOCATE(a)
		DEALLOCATE(Ts)
		DEALLOCATE(w)

		!try to get the particular solution
		DO i=1,n
			A0(i,1:m)= MATMUL((v(NumDs,1:m,1:m)), u(i,1:m)*Uw(1:m))
		END DO
		
		DO i=1,m
			DO j=1,m
				DO q=1, n
					Jk(NumDs,i,j)=Jk(NumDs,i,j) + (DifX(i,q)-B(i,q))*A0(q,j)
				END DO
			END DO
		END DO

		! It's time to release the memory 
		DEALLOCATE(B)
		DEALLOCATE(DifX)  
		DEALLOCATE(Uw)
		DEALLOCATE(u)
		DEALLOCATE(A0)	
	END DO
	CLOSE(7)

	!Here we need to decide which weight system should be used
	Weight= NumSample/NumSapAll
	!Weight=Wdata         !This weight is got from the data file

	!The thrid part: iterated to generate the common solution
	!Computing the initial value
	ALLOCATE(J0(NumDs,m,m))
	ALLOCATE(Y(NumDs,m,m))
	ALLOCATE(OptGNW(m,m))
	Y=0.0
	OptGNW=0.0
	J0=JK
	GNW=0.0

	!computing the object function	
	ObjFun=0.0
	DO k=1,NumDs
		DO i=1,m
			DO j=1,m
				ObjFun = ObjFun+ (Weight(k)*ABS(GNW(i,j) - Jk(k,i,j)) +lamda*ABS(GNW(i,j))) /(m**2*NumDS)  !negelect the Y(k)
			END DO
		END DO
	END DO

	PreObjFun=100000.0
	IteNum=0

	DO WHILE ((ABS(PreObjFun-ObjFun) > epsilon) .AND. (IteNum<IteMax))
		IteNum= IteNum+1
		WRITE(8,*) '***********************************'
		WRITE(8, '("The Iteration number is ", I4)')  IteNum		
		WRITE(8,*) '***********************************'
		PreObjFun=ObjFun
		OptGNW=GNW

		!Step 1: Fixing the J, using L1 regression to find proper solution of Yk
		DO k=1,NumDs
			IF (L(k)>0) THEN
				n=NumSample(k)-1

				!Dynamic allocate the memory space for the arrays
				DO dc=1, m
					ALLOCATE(L1S(m))
					ALLOCATE(L1A(m+2, L(k)+2))
					ALLOCATE(L1B(m))
					ALLOCATE(L1E(m))
					ALLOCATE(L1X(L(k)))
					L1S=0.0
					L1A=0.0
					L1B=0.0
					L1E=0.0
					L1X=0.0

					! Solving the L1 regression problem respectively	
					DO j=1,m
						L1B(j)=  GNW(dc,j)-J0(k,dc,j)  
					END DO
					
					DO i=1,m
						DO j=1, L(k)
							L1A(i,j)=v(k,i,j)
						END DO
					END DO

					TOLER=10**(-3*2/3)
					!write(*, '("Begin to solve the decomposited subproblem  ", I4)')  dc
					CALL L1(m,L(k),m+2,L(k)+2,L1A,L1B,TOLER,L1X,L1E,L1S)
					! A(M+1,N+1)  THE MINIMUM SUM OF THE ABSOLUTE VALUES OF THE RESIDUALS.
					! A(M+1,N+2)  THE RANK OF THE MATRIX OF COEFFICIENTS.
					! A(M+2,N+1)  EXIT CODE WITH VALUES.
					!             0 - OPTIMAL SOLUTION WHICH IS PROBABLY NON-
					!                 UNIQUE (SEE DESCRIPTION).
					!             1 - UNIQUE OPTIMAL SOLUTION.
					!             2 - CALCULATIONS TERMINATED PREMATURELY DUE TO
					!                 ROUNDING ERRORS.
					! A(M+2,N+2)  NUMBER OF SIMPLEX ITERATIONS PERFORMED.
					!write(*, '("Return values of the L1 solver ", f12.6,f12.6,f12.6,f12.6)')  L1A(m+1,L(k)+1),L1A(m+1,L(k)+2),L1A(m+2,L(k)+1),L1A(m+2,L(k)+2)
		
					DEALLOCATE(L1S)
					DEALLOCATE(L1A)
					DEALLOCATE(L1B)
					DEALLOCATE(L1E)
					
					!Give the solution to the Y matrix
					DO j=1,m
						IF (j> L(k)) THEN
							Y(k, dc, j)=0
						ELSE
							Y(k, dc, j)=L1X(j)
						END IF
					END DO
					DEALLOCATE(L1X)
				END DO
			END IF
		END DO


		! Step 2: Fixing the Y, solve the value of J	
		DO k=1,NumDs
			DO i=1, m
				DO j=1,m
					DO q=1,m
						J0(k,i,j)=J0(k,i,j)+Y(k,i,q)*v(k,j,q)
					END DO
				END DO
			END DO
		END DO
		
		IF (NumDs==1) THEN
			GNW=J0(NumDs,1:m,1:m)
			EXIT
		END IF

		!computing the object function	
		objFun=0
		DO q=1,NumDs
			DO i=1,m
				DO j=1,m
					ObjFun = ObjFun+ (Weight(q)*ABS(GNW(i,j) - J0(q,i,j)) +lamda*ABS(GNW(i,j))) /(m**2*NumDS)
				END DO
			END DO
		END DO

		!Trying to give the exact solution
		ALLOCATE(SortRoot(NumDs+1))
		ALLOCATE(Omega(NumDs+1))
		SortRoot=0.0
		Omega=0.0

		DO i=1,m
			DO j=1,m
				DO k=1,NumDs
					SortRoot(k)= J0(k,i,j)
					Omega(k)=Weight(k)
				END DO

				SortRoot(NumDs+1)=0.0
				Omega(NumDs+1)=NumDs*lamda

				!Sorting the root in increasing order
				DO k=1,NumDs+1
					Location=MAXLOC(SortRoot(1:NumDs-k+2)) 
					temp=SortRoot(NumDs-k+2)
					SortRoot(NumDs-k+2) = SortRoot(Location(1))
					SortRoot(Location(1))=temp

					temp=Omega(NumDs-k+2)
					Omega(NumDs-k+2) = Omega(Location(1))
					Omega(Location(1))=temp
				END DO

				CALL SumAbsMin(NumDs+1, SubMin, SubObj, Omega(1:NumDs+1), SortRoot(1:NumDs+1))

				GNW(i,j)=subMin   
			END DO
		END DO
		DEALLOCATE(SortRoot)
		DEALLOCATE(Omega)

		!computing the object function	
		objFun=0
		DO k=1,NumDs
			DO i=1,m
				DO j=1,m
					ObjFun = ObjFun+ (Weight(k)*ABS(GNW(i,j) - J0(k,i,j)) +lamda*ABS(GNW(i,j))) /(m**2*NumDS)
				END DO
			END DO
		END DO
		WRITE(8, '("The values of object funciton is ", f20.10 )')  ObjFun
	END DO

	OptGNW=GNW
	!print the final results
	WRITE(8,*) 'The final strength and influence of genes matrix'
	DO i=1,m
		WRITE(8,*) OptGNW(i,1:m)
	END DO

	! make the statistical report
	Nzero=0
	DO i=1,m
		DO j=1,m
			IF (abs(OptGNW(i,j)) < Threshold) THEN
				Nzero=Nzero+1
			END IF
		END DO
	END DO

	WRITE(8, '("The number of the non zero element ", I8, f12.6)')  m*m-Nzero, 1-real(Nzero)/m**2

	ALLOCATE(Var(m,m))
	ALLOCATE(Sigma(m,m)) 
	Var=0.0
	Sigma=0.0
	AveSigma=0
	AveV=0
	DO i=1,m
		DO j=1,m
			DO k=1,NumDs
				Var(i,j) = Var(i,j) + (Weight(k)*(GNW(i,j) - J0(k,i,j))**2) /(NumDS)
			END DO
			Sigma(i,j) = sqrt(Var(i,j))
			AveSigma= AveSigma +Sigma(i,j)
			AveV=AveV + Var(i,j)
		END DO
	END DO
	AveSigma= AveSigma/(m**2)
	AveV=AveV/(m**2)
	WRITE(8, '("The Confidence Evaluation ", f12.6 , f12.6)') AveV, AveSigma
	CLOSE(8)

	!Writing the inferred network to file
	CALL Write_Network(inputFile,OptGNW, m)
	CALL Write_DotFile(inputFile,GeneName, OptGNW, m, Threshold)

	DEALLOCATE(J0)
	DEALLOCATE(GNW)
	DEALLOCATE(OptGNW)
	DEALLOCATE(Y)
	DEALLOCATE(v)
	DEALLOCATE(JK)
	DEALLOCATE(Var)
	DEALLOCATE(Sigma)
	DEALLOCATE(GeneName)
END SUBROUTINE mainprocess
