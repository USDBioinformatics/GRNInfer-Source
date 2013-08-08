	SUBROUTINE SumAbsMin(Dim, Submin, Subobj, Weight, SortRoot)
	! This subroutine aims to solve the subproblem to get J(ij)
	! The input is the sorted root
	IMPLICIT NONE	
	INTEGER(4),  INTENT(IN) :: Dim
	REAL(8), DIMENSION(Dim), INTENT(IN) :: Weight, SortRoot
	REAL(8), INTENT(OUT) :: SubMin, SubObj
	
	INTEGER(4) :: i,j,k,q
	REAL(8) :: tempVal, tempWei, SubMinTemp, SubObjTemp,SumOmega 
				! Give the inital value as the all negative case
				tempVal=0.0
				tempWei=0.0
				do k=1, Dim
					tempVal= tempVal + weight(k)*SortRoot(k)
					tempWei= tempWei + weight(k)
				end do

				SubMin=tempVal/tempWei
				SubObj=0.0
				if (SubMin>SortRoot(1)) then
					SubMin= SortRoot(1)
					SubObj=  tempVal- tempWei*SortRoot(1)
				end if
				
				! Iterate through other segment to find the real solution
				do k=1, Dim-1
					tempVal=0.0
					tempWei=0.0
					do q=1, Dim
						if (q>=k+1) then
							tempVal= tempVal + weight(q)*SortRoot(q)
							tempWei= tempWei - weight(q)
						else
							tempVal= tempVal - weight(q)*SortRoot(q)
							tempWei= tempWei + weight(q)
						end if			
					end do
					
					SubMinTemp=-tempVal/tempWei    				
					if ((SubMinTemp > SortRoot(k)) .AND. (SubMinTemp <= SortRoot(k+1))) then
						SubObjTemp=  0.0
					else
						if (tempWei < 0.0) then   
							SubMinTemp= SortRoot(k+1)
							SubObjTemp=  tempVal + tempWei*SubMinTemp     
						else       
							SubMinTemp= SortRoot(k)
							SubObjTemp=  tempVal + tempWei*SubMinTemp    
						end if			
					end if

					if (SubObjTemp < SubObj) then
						SubMin= SubMinTemp
						SubObj= SubObjTemp
					end if

				end do	

				!Last Consider the case that all positive				
				tempVal=0.0
				tempWei=0.0
				do k=1, Dim
					tempVal= tempVal - weight(k)*SortRoot(k)
					tempWei= tempWei + weight(k)
				end do
				SubMinTemp= -tempVal/tempWei
				SubObjTemp= 0.0

				if (SubMinTemp < SortRoot(Dim)) then
					SubMinTemp = SortRoot(Dim)
					SubObjTemp = tempVal + tempWei*SortRoot(Dim)
				end if

				if (SubObjTemp < SubObj) then
					SubMin= SubMinTemp
					SubObj= SubObjTemp
				end if

END SUBROUTINE SumAbsMin