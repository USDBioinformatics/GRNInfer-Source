	SUBROUTINE Write_Network(inputFile,Network, Dim)

	IMPLICIT NONE	
	INTEGER(4),  INTENT(IN) :: Dim
	REAL(8), DIMENSION(Dim,Dim), INTENT(IN) :: Network
	INTEGER(4) :: i,j
	CHARACTER(LEN=*) :: inputFile

	OPEN (UNIT=3, FILE=inputFile//'.dat', STATUS='REPLACE')
	
	do i=1,Dim
		do j=1, Dim
			write(3,'(f25.12 )') Network(i,j)
		end do
		write (3,'(:)')
	end do
	close(3)
	END SUBROUTINE Write_network


	SUBROUTINE Write_DotFile(inputFile,GeneName,Network,Dim,Threshold)

	IMPLICIT NONE	
	INTEGER(4),  INTENT(IN) :: Dim
	CHARACTER*9, DIMENSION(Dim), INTENT(IN) :: GeneName
	REAL(8), DIMENSION(Dim,Dim), INTENT(IN) :: Network
	double precision, INTENT(IN) :: Threshold
	INTEGER(4) :: i,j
	CHARACTER(LEN=*) :: inputFile;

	OPEN (UNIT=5, FILE=inputFile//'.dot', STATUS='REPLACE')

	write(5,*) 'digraph "Gene regulatory network" {'
	write(5,*) '		graph                            '
	write(5,*) '		[                                '
	write(5,*) '	center="true"                        '
	write(5,*) '			overlap="false"              '
	write(5,*) '			Damping=0.999                '
	write(5,*) '			fontname="Helvetica"         '
	write(5,*) '			maxiter=1000000              '
	write(5,*) '			splines="true"               '
	write(5,*) '			sep=0.8                      '
	write(5,*) '			epsilon=0.0000001            '
	write(5,*) '			ratio="auto"                 '
	write(5,*) '		]                                '
	write(5,*) '		node                             '
	write(5,*) '		[                                '
	write(5,*) '			fontsize=20                  '
	write(5,*) '			fontname="Helvetica-bold"    '
	write(5,*) '			shape="circle"               '
	write(5,*) '			style="bold"                 '
	write(5,*) '		]                                '
	write(5,*) '		edge                             '
	write(5,*) '		[                                '
	write(5,*) '			fontsize=18                   '
	write(5,*) '			fontname="Helvetica"         '
	write(5,*) '			arrowsize=2.0                 ' 
	write(5,*) '			len=2.5                      '  
	write(5,*) '		]                                '

	do i=1,Dim
		do j=1, Dim
			if (Network(i,j) > Threshold) then        !The activation is denoted by blue
				write(5,'( A9 "  ->  ", A9  )') , GeneName(j), GeneName(i)
				write(5,*) '		[color="red", arrowhead="normal"]'
			elseif (Network(i,j) < -Threshold) then   !The repressison is denoted by red
				write(5,'( A9 "  ->  ", A9  )') , GeneName(j), GeneName(i)
				write(5,*) '		[color="blue", arrowhead="tee"] '
			end if
		end do
	end do
	write(5,*) '	} '
	close(5)

	END SUBROUTINE Write_DotFile
