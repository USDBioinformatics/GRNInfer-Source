MODULE large_array

	!Dynamic allocating the memory then can release them when necessayr
	double precision, allocatable :: a(:,:), B(:,:),DifX(:,:)  !DIMENSION(ND,MP,NP) :: a, B, DifX 
	double precision, allocatable :: Ts(:)                     !DIMENSION(ND,NP) :: Ts
	double precision, allocatable :: Uw(:), w(:)               !DIMENSION(ND,MP) :: w, Uw
	double precision, allocatable :: u(:,:)                    !DIMENSION(ND,NP,MP) :: u
	double precision, allocatable :: A0(:,:)                   !DIMENSION(NP,MP) :: A0


	double precision, allocatable::   J0(:,:,:), Y(:,:,:), JK(:,:,:), v(:,:,:)    !DIMENSION(ND,MP,MP) :: v, J0, Y

	double precision, allocatable:: GNW(:,:), OptGNW(:,:), Var(:,:), Sigma(:,:)   !DIMENSION(MP,MP) :: GNW,Jk
	double precision, allocatable  ::   sortRoot(:), Omega(:)

	!Dynamic Allocating the memory for these arrays
	INTEGER, allocatable :: L1S(:)                                !DIMENSION(MP) :: L1S
	double precision, allocatable  ::   L1A(:,:)                  !DIMENSION(MP+2,MP+2) :: L1A
	double precision, allocatable  ::   L1B(:), L1E(:), L1X(:)    !DIMENSION(MP) :: L1B,L1E, L1X
	                                                              

	!To save the gene name 
	CHARACTER*9, allocatable :: GeneName(:) 

END MODULE large_array
