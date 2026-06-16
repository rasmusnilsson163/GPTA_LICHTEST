! ! Code written by Rasmus Nilsson and implemented into GPTA written by
! ! copyrightholder Paolo Raiteri. The implemented theory is based on
! ! the LICH-TEST method presented in https://pubs.acs.org/doi/10.1021/acs.jpcb.1c01926 

! ! This program is free software; you can redistribute it and/or modify it 
! ! under the terms of the GNU General Public License as published by the 
! ! Free Software Foundation; either version 3 of the License, or 
! ! (at your option) any later version.
! !  
! ! Redistribution and use in source and binary forms, with or without 
! ! modification, are permitted provided that the following conditions are met:
! ! 
! ! * Redistributions of source code must retain the above copyright notice, 
! !   this list of conditions and the following disclaimer.
! ! * Redistributions in binary form must reproduce the above copyright notice, 
! !   this list of conditions and the following disclaimer in the documentation 
! !   and/or other materials provided with the distribution.
! ! * Neither the name of the <ORGANIZATION> nor the names of its contributors 
! !   may be used to endorse or promote products derived from this software 
! !   without specific prior written permission.
! ! 
! ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! ! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! ! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! ! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! ! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! ! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! ! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! ! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! ! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! ! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! !
module moduleLichTest
  use moduleVariables
  use moduleFiles
  use moduleSystem
  use moduleMessages
  use moduleStrings

  implicit none

  public :: computeLichTest, computeLichTestHelp
  private
  
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  
  integer, pointer :: tallyExecutions
  real(8), pointer :: rcut
  real(8), pointer :: min_score
  integer, pointer, dimension(:) :: labels

contains

  subroutine computeLichTestHelp()
    use moduleMessages
    implicit none
    call message(0,"This action uses LICH-TEST to classify ice structures")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --lichtest +s Ow,o +out lichtest.out")
    call message(0,"  gpta.x --i coord.pdb --lichtest +s mW --o labeled_ice.pdb")
    call message(0,"  gpta.x --i coord.pdb --lichtest +s mW +minscore 0.6 --o labeled_ice.pdb")
    call message(0,"  gpta.x --i coord.pdb --lichtest +s o,oh,ohs +rcut 2.5 --o labeled_ice.pdb")
  end subroutine computeLichTestHelp

  subroutine computeLichTest(a)
    use moduleMessages
    implicit none
    type(actionTypeDef), target :: a

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      if (firstAction) then
        call dumpScreenInfo()

        if (resetFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if

        ! select one group of atoms
        call selectAtoms(1,actionCommand,a)
        ! create a list of the atoms' indices for each group
        call createSelectionList(a,1)

        ! Check groups have atoms
        if (count(a % isSelected(:,1)) == 0) call message(-1,"--lichtest - no atoms selected")

        ! Throw a warning for unused flags
        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else

        if (a % updateAtomsSelection) then 
          call selectAtoms(1,actionCommand,a)
          call createSelectionList(a,1)
        end if

      end if

      allocate(labels(frame % natoms))

      call computeAction(a)

    end if

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) call finaliseAction(a)

  end subroutine computeLichTest

  subroutine initialiseAction(a)
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a

    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile
    rcut                 => a % doubleVariables(1)
    min_score            => a % doubleVariables(2)

    a % actionInitialisation = .false.
    a % cutoffNeighboursList = 3.5d0

    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .true.
    a % requiresNeighboursListDouble = .true.

    call assignFlagValue(actionCommand,"+rcut",rcut,3.5d0)
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'lich_test.out')
    call assignFlagValue(actionCommand,"+minscore",min_score,0.5d0)

    call initialiseFile(outputFile,outputFile % fname)                
    write(outputFile % funit,"(a)") "# Frame Number |       Liquid |        Cubic |    Hexagonal |        Mixed |    Cubic-if. |Hexagonal-if. | Clath. hydr. |  Interfacial | Clathrate-if."
  
    !a % cutoffNeighboursList = rcut
    
    tallyExecutions = 0 

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    use moduleMessages
    implicit none
    call message(0,"Settings for LICH-TEST")
    call message(0,"...Minimal Score   ",r=min_score)
    call message(0,"...Cutoff distance ",r=rcut)
  
  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleVariables
    use moduleSystem
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a

    ! Matrices
    integer :: P(6,3,3)
    real(8) :: Ts(3,3), Te(3,3), Tx(3,3)
    real(8) :: U(4,3), V(4,3)
    real(8), allocatable, dimension(:,:,:) :: neighbour_vectors
    integer, allocatable, dimension(:,:) :: neighbour_indices
    integer, allocatable, dimension(:,:) :: neigh_labels

    ! Arrays
    real(8) :: dij(3)
    integer, allocatable, dimension(:) :: neighbour_number

    ! Other 
    integer, allocatable, dimension(:,:) :: S, E
    real(8), allocatable, dimension(:,:) :: S_real, E_real
    integer :: nsel
    
    ! Helpers
    integer :: iatm, ineigh, jatm
    integer :: i, j, k 
    integer :: column_pos, discon_flag 
    integer :: nU, nV 
    real(8) :: s_scr, e_scr
    real(8) :: dist1, dist2, dist3, dist4, dist_temp

    nsel = count( a % isSelected(:,1))
    ! Allocate arrays
    allocate(S(frame % natoms,4),E(frame % natoms,4) , source=0)
    allocate(S_real(frame % natoms,4),E_real(frame % natoms,4) , source=0.0d0)
    allocate(neighbour_vectors(frame % natoms,3,4)) ! All atoms with 4 neighbours with 3 coordinates each
    allocate(neighbour_indices(frame % natoms,4))   ! All atoms has 4 neigbour with indices
    allocate(neighbour_number(frame % natoms))      ! How many close neighbours an atom have, integer from 0 to 4
    
    neighbour_vectors = 0.0d0
    neighbour_indices = -1
    neighbour_number = 0
    S = 0
    E = 0
    S_real = 0.0d0
    E_real = 0.0d0
        
    ! Define template matrices
    Ts(1,:) = (/ -1.0d0,  0.5d0,  0.5d0 /) 
    Ts(2,:) = (/  0.5d0, -1.0d0,  0.5d0 /) 
    Ts(3,:) = (/  0.5d0,  0.5d0, -1.0d0 /)

    Te(1,:) = (/  0.5d0, -0.5d0, -0.5d0 /)
    Te(2,:) = (/ -0.5d0,  0.5d0, -0.5d0 /)
    Te(3,:) = (/ -0.5d0, -0.5d0,  0.5d0 /)
    
    P(1,:,:) = reshape((/ 0, 0, 1, 0, 1, 0, 1, 0, 0 /), (/3,3/))
    P(2,:,:) = reshape((/ 0, 0, 1, 1, 0, 0, 0, 1, 0 /), (/3,3/))
    P(3,:,:) = reshape((/ 0, 1, 0, 0, 0, 1, 1, 0, 0 /), (/3,3/))
    P(4,:,:) = reshape((/ 0, 1, 0, 1, 0, 0, 0, 0, 1 /), (/3,3/))
    P(5,:,:) = reshape((/ 1, 0, 0, 0, 0, 1, 0, 1, 0 /), (/3,3/))
    P(6,:,:) = reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /), (/3,3/))    
 

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,iatm,ineigh,jatm,dij,column_pos,dist1,dist2,dist3,dist4,dist_temp)
!$OMP DO
      ! Create neighbour_vectors
      do i = 1, nsel
        iatm = a % idxSelection(i,1)
        dist1 = 3.5d0; dist2 = 3.5d0; dist3 = 3.5d0; dist4 = 3.5d0
        ! Loop over iatom neighbours, first find 4 closest to create V matrix
        do ineigh = 1, nneigh(iatm)

          if (rneigh(ineigh,iatm) > rcut) cycle

          ! Find index and vector of neigbour atom
          jatm = lneigh(ineigh,iatm)
          if ( .not. a % isSelected(jatm,1) ) cycle
         ! write(*,*) frame % pos(:,jatm) ,',', frame % pos(:,iatm)
          dij = frame % pos(:,jatm) - frame % pos(:,iatm)
          dist_temp = computeDistanceSquaredPBC(dij)
          dist_temp = sqrt(dist_temp)

          ! Check if distance in at least 4th closest
          if (dist_temp < dist1) then
            column_pos = 1
            dist4 = dist3; dist3 = dist2; dist2 = dist1; dist1 = dist_temp
          else if (dist_temp < dist2) then
            column_pos = 2
            dist4 = dist3; dist3 = dist2; dist2 = dist_temp
          else if (dist_temp < dist3) then
            column_pos = 3
            dist4 = dist3; dist3 = dist_temp
          else if (dist_temp < dist4) then
            column_pos = 4
            dist4 = dist_temp
          else
            cycle
          end if
          
          ! Valid neighbour, add 1 if not already 4
          if (neighbour_number(iatm) < 4) then
            neighbour_number(iatm) = neighbour_number(iatm) + 1
          end if

          ! Normalize vector
          dij = dij / sqrt(sum(dij*dij))
          
          ! Shift columns to the right, thus the closest atom is always in 
          ! the left-most column
          do j = 3, column_pos, -1
            neighbour_vectors(iatm,:,j+1) = neighbour_vectors(iatm,:,j)
            neighbour_indices(iatm,j+1) = neighbour_indices(iatm,j)
          end do
          neighbour_vectors(iatm,:,column_pos) = dij
          neighbour_indices(iatm,column_pos) = jatm

        end do
      end do
      ! All neigbour matices made now
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (iatm,j,i,jatm,Tx,nU,nV,discon_flag,s_scr,e_scr,k)
!$OMP DO
      do i = 1, nsel
        iatm = a % idxSelection(i,1)
        do j = 1,4
          jatm = neighbour_indices(iatm,j)
          if (jatm > iatm) then
            Tx = 0
            ! In neigbour vectors, 1st which atom, 2nd atom coords, 3rd neigh atom
            call UtV_func(neighbour_vectors(iatm,:,:),neighbour_vectors(jatm,:,:),Tx,discon_flag) 
            nU = neighbour_number(iatm) - 1
            nV = neighbour_number(jatm) - 1
        
            if (discon_flag == 0) then
              call score_func(Tx,Ts,Te,nU,nV,P,s_scr,e_scr) 
            else 
              neighbour_indices(iatm,j) = 0
              neighbour_number(iatm) = neighbour_number(iatm) - 1
              s_scr = 0.0d0
              e_scr = 0.0d0
            end if
            
            ! The following code coresponds to E(E<S)=0; S(S<E)=0; E(E<score_min)=0; S(S<score_min)=0; 
            E_real(iatm,j) = e_scr
            S_real(iatm,j) = s_scr

            do k = 1,4
              if (neighbour_indices(jatm,k) .eq. iatm) then
                S_real(jatm,k) = s_scr
                E_real(jatm,k) = e_scr
              end if
            end do

          end if
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
    deallocate(neighbour_vectors)
    
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,iatm,j)
!$OMP DO
    ! The following code coresponds to E(E<S)=0; S(S<E)=0; E(E<score_min)=0; S(S<score_min)=0; 
    do i = 1,nsel
      iatm = a % idxSelection(i,1)
      do j = 1,4
        if (E_real(iatm,j) .lt. S_real(iatm,j)) E_real(iatm,j) = 0.0d0
        if (S_real(iatm,j) .lt. E_real(iatm,j)) S_real(iatm,j) = 0.0d0
        if (E_real(iatm,j) .lt. min_score)      E_real(iatm,j) = 0.0d0
        if (S_real(iatm,j) .lt. min_score)      S_real(iatm,j) = 0.0d0 
      
        ! Round values
        E(iatm,j) = nint(E_real(iatm,j))
        S(iatm,j) = nint(S_real(iatm,j))
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! Now to label the ice
    allocate(neigh_labels(frame % natoms,4))
    neigh_labels = 0
    labels = 0

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,j,iatm)
!$OMP DO
      ! First three bullet points from the article
      do i = 1, nsel
        iatm = a % idxSelection(i,1)
        if ( (sum(S(iatm,:)) .eq. 4) .and. (sum(E(iatm,:)) .eq. 0) ) labels(iatm) = 1 ! Cubic
        if ( (sum(S(iatm,:)) .eq. 3) .and. (sum(E(iatm,:)) .eq. 1) ) labels(iatm) = 2 ! Hexagonal
        if (  sum(E(iatm,:)) .ge. 4) labels(iatm) = 6 ! Clathrate hydrate
      end do
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,j,iatm)
!$OMP DO
      ! Interfacial labeling
      do i = 1, nsel
        iatm = a % idxSelection(i,1)
      !  print *, sum(S(iatm,:)) + sum(E(iatm,:))
        if ( (sum(S(iatm,:)) + sum(E(iatm,:)) .eq. 0) .or. (labels(iatm) .gt. 0) ) cycle
        ! Find possible label on neigbours
        do j = 1,  neighbour_number(iatm)
          neigh_labels(iatm,j) = labels(neighbour_indices(iatm,j))
        end do
        ! The remaining 5 points in the article are tested here
        if ( (any( neigh_labels(iatm,:) .eq. 1 ) .and. any( neigh_labels(iatm,:) .eq. 2)) .or. &
          &  (any( neigh_labels(iatm,:) .eq. 1 ) .and. any( neigh_labels(iatm,:) .eq. 6)) .or. &
          &  (any( neigh_labels(iatm,:) .eq. 2 ) .and. any( neigh_labels(iatm,:) .eq. 6)) ) then
          labels(iatm) = 3 ! Mixed
        else if ( any( neigh_labels(iatm,:) .eq. 1)) then
          labels(iatm) = 4 ! CI, cubic-interfacial
        else if ( any( neigh_labels(iatm,:) .eq. 2)) then
          labels(iatm) = 5 ! HI, hexagoanl-interfacial
        else if ( any( neigh_labels(iatm,:) .eq. 6)) then
          labels(iatm) = 8 ! CHI, clathrate-interfacial 
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,j,iatm)
!$OMP DO
      ! Interfacial labeling
      do i = 1, nsel
        iatm = a % idxSelection(i,1)
      !  print *, sum(S(iatm,:)) + sum(E(iatm,:))
        if ( (sum(S(iatm,:)) + sum(E(iatm,:)) .eq. 0) .or. (labels(iatm) .gt. 0) ) cycle
        ! Find possible label on neigbours
        do j = 1,  neighbour_number(iatm)
          neigh_labels(iatm,j) = labels(neighbour_indices(iatm,j))
        end do
        
        if ( (any( neigh_labels(iatm,:) .eq. 3 )) .or. &
          &  (any( neigh_labels(iatm,:) .eq. 4 )) .or. &
          &  (any( neigh_labels(iatm,:) .eq. 5 )) .or. &
          &  (any( neigh_labels(iatm,:) .eq. 8 ))) then
          labels(iatm) = 7 ! I, other interfacial ice 
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL


    deallocate(S,E,S_real,E_real,neighbour_number,neigh_labels,neighbour_indices)
    
    ! Set new labels to atoms in frame
    call set_labels(a)
    call write_label_counts(a)

  end subroutine computeAction

  subroutine score_func(Tx,Ts,Te,nU,nV,P,s_scr,e_scr)
    implicit none
    real(8), intent(inout) :: s_scr, e_scr
    real(8), intent(in) :: Tx(:,:), Ts(:,:), Te(:,:)
    integer, intent(in) :: P(:,:,:)
    integer, intent(in) :: nU, nV
    ! Param
    real(8) :: lmb = 0.15d0
    ! Helpers
    integer :: i,j,t
    real(8) :: tmp_s_scr, tmp_e_scr

    tmp_s_scr = 0.0d0
    tmp_e_scr = 0.0d0

    if ((nU .gt. 0) .and. (nV .gt. 0)) then ! Else the values are 0
    ! t = 1
    t = 1
    do i = 1,nV
      do j = 1,3                               ! D^2(i,j)             P_t(i,j)
        tmp_s_scr = tmp_s_scr + ( sum((Tx(1:nU,i) - Ts(1:nU,j))**2) * P(t,i,j) )
        tmp_e_scr = tmp_e_scr + ( sum((Tx(1:nU,i) - Te(1:nU,j))**2) * P(t,i,j) )
      end do
    end do
    s_scr = tmp_s_scr
    e_scr = tmp_e_scr
    
    ! t = 2 to 6
    do t = 2,6
      tmp_s_scr = 0.0d0
      tmp_e_scr = 0.0d0
      do i = 1,nV
        do j = 1,3
          tmp_s_scr = tmp_s_scr + ( sum((Tx(1:nU,i) - Ts(1:nU,j))**2) * P(t,i,j) )
          tmp_e_scr = tmp_e_scr + ( sum((Tx(1:nU,i) - Te(1:nU,j))**2) * P(t,i,j) )
        end do
      end do
      if (tmp_s_scr < s_scr) s_scr = tmp_s_scr
      if (tmp_e_scr < e_scr) e_scr = tmp_e_scr
    end do 

    s_scr = exp(-1.0d0*s_scr / (lmb*nU*nV)) ! = S_s
    e_scr = exp(-1.0d0*e_scr / (lmb*nU*nV)) ! = S_e
    
    else
    
    s_scr = 0.0d0
    e_scr = 0.0d0

    end if
    
  end subroutine score_func

  subroutine UtV_func(U,V,Tx_out,disconnection_flag)
    implicit none
    real(8), intent(in) :: U(:,:), V(:,:) ! Both 3 x 4 matrices
    real(8), intent(inout) :: Tx_out(:,:)
    integer, intent(inout) :: disconnection_flag

    ! Helpers
    integer :: i, j, k, l, ii, jj
    real(8) :: Tx(4,4)
    
    ! NB! The columns in the U and V matrixes are the atoms,
    ! the rows are the positions of them,
    ! i.e, V(a1,a2,a3,a4) = V(a1(x,y,z),a2(x,y,z),...,a4(x,y,z)) 

    ! ! U and V matrixes columns need to be normalized. This should be already done
    ! Compute U^T * V
    Tx = matmul(transpose(U),V)
   
    ! Find indices where abs(Tx + 1) < 1e-7
    ii = 0
    jj = 0
    do i = 1, 4 ! size(Tx, 1)
        do j = 1, 4 ! size(Tx, 2)
            if (abs(Tx(i, j) + 1.0) < 1.0e-7) then
                ii = i
                jj = j
                exit
            end if
        end do
    end do

    ! Check if there is no central O-O bond
    if (ii == 0) then
      disconnection_flag = 1
    else
      ! Remove column and row where -1 appears
      k = 1
      do i = 1, 4
          if (i /= ii) then
              l = 1
              do j = 1, 4
                  if (j /= jj) then
                      Tx_out(k, l) = Tx(i, j)
                      l = l + 1
                  end if
              end do
              k = k + 1
          end if
      end do
      disconnection_flag = 0
    end if
  end subroutine UtV_func

  subroutine write_label_counts(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: i, iatm
    integer :: label_counts(9)

    label_counts = 0
 
    ! Write out file 
    do i = 1, count( a % isSelected(:,1))
      iatm = a % idxSelection(i,1)
      label_counts(labels(iatm)+1) = label_counts(labels(iatm)+1) + 1
    end do

    write(outputFile % funit,"(10(I15))") frame % nframe, label_counts(1), label_counts(2), label_counts(3), label_counts(4), label_counts(5), label_counts(6), label_counts(7), label_counts(8), label_counts(9)
  end subroutine write_label_counts

  subroutine set_labels(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: i, iatm
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,iatm)
!$OMP DO
      do i = 1, count( a % isSelected(:,1))
        iatm = a % idxSelection(i,1)
        if      (labels(iatm) .eq. 0) then
          frame % lab(iatm) = "OL  "  ! Liquid
        else if (labels(iatm) .eq. 1) then
          frame % lab(iatm) = "OC  "  ! Cubic
        else if (labels(iatm) .eq. 2) then
          frame % lab(iatm) = "OH  "  ! Hexagonal
        else if (labels(iatm) .eq. 3) then
          frame % lab(iatm) = "OM  "  ! Mixed interfacial
        else if (labels(iatm) .eq. 4) then
          frame % lab(iatm) = "OIC "  ! Cubic interfacial
        else if (labels(iatm) .eq. 5) then
          frame % lab(iatm) = "OIH "  ! Hexagonal interfacial
        else if (labels(iatm) .eq. 6) then
          frame % lab(iatm) = "OCH "  ! Clathrate hydrate
        else if (labels(iatm) .eq. 7) then
          frame % lab(iatm) = "OI  "  ! Other interfacial
        else if (labels(iatm) .eq. 8) then
          frame % lab(iatm) = "OICH"  ! Clathrade interfacial
        end if 
      end do
    
!$OMP END DO
!$OMP END PARALLEL
  end subroutine set_labels


  subroutine finaliseAction(a)
    implicit none
    type(actionTypeDef), target :: a
   
    close(outputFile % funit)
    
  end subroutine finaliseAction


end module moduleLichTest