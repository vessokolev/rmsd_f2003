subroutine s_rmsd_superimpose(num_atoms,do_count,x1,x2,rmsd)
!
! This subroutine compares two sets of points (atoms) of the same size (!!!),
! computes RMSD, and superimposes the two coordinate sets - by translating and
! rotating the set "x2". Do back up the array "x2" to produce a copy of the
! previous values, in case they need to be restored back later.
!
! This subroutine is part of the implementation of Quaternion Characteristic
! Polynomial (QCP) method, based on the code published by Naoto Hori at:
!
! https://github.com/naotohori/fQCP
!
! under BSD license. His code, in turn, is based on the C code written by Pu Liu
! and Douglas L. Theobald. The code bellow is improved and partially revised by
! Veselin Kolev <vesso.kolev@gmail.com>.
!
! To cite QCP:
!
! 1. Douglas L. Theobald (2005) "Rapid calculation of RMSD using a quaternion-based
!    characteristic polynomial." Acta Crystallographica A 61(4):478-480.
!
! 2. Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2010) "Fast
!    determination of the optimal rotational matrix for do_counted
!    superpositions." Journal of Computational Chemistry 31(7):1561-1563
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! License: BSD
! Version: 2018121200
!
use iso_c_binding,only:C_INT,C_FLOAT
!
implicit none 
!
! Interface variables:
!
integer(C_INT),intent(in)   :: num_atoms
logical,intent(in)          :: do_count(num_atoms)
real(C_FLOAT),intent(in)    :: x1(3,num_atoms)
real(C_FLOAT),intent(inout) :: x2(3,num_atoms)
real(C_FLOAT),intent(out)   :: rmsd
!
! Local variables:
!
integer(C_INT)              :: i
integer(C_INT)              :: j
integer(C_INT)              :: ireturn
integer(C_INT)              :: ncount
integer(C_INT),parameter    :: t(6,3)=reshape(source=(/1,1,2,2,3,3,4,1,5,2,6,3,7,1,8,2,9,3/),shape=(/6,3/))
real(C_FLOAT)               :: A(9)
real(C_FLOAT)               :: E0
real(C_FLOAT)               :: rot(9)
real(C_FLOAT)               :: x1_bak(3,num_atoms)
real(C_FLOAT)               :: x2_bak(3,num_atoms)
real(C_FLOAT)               :: center1(3)
real(C_FLOAT)               :: center2(3)
real(C_FLOAT)               :: x_tmp(3)

! Create a copy of the input coordinates:
!
x1_bak=x1(:,:)
x2_bak=x2(:,:)
!
! Center the atom coordinates for the first molecule:
!
call s_rmsd_center_coords(num_atoms,do_count,x1_bak,ncount)
!
! Center the atom coordinates for the second molecule:
!
call s_rmsd_center_coords(num_atoms,do_count,x2_bak,ncount)
!
! Compute the inner product of the coordinates of the molecules:
!
call s_rmsd_inner_product(num_atoms,do_count,x1_bak,x2_bak,A,E0)
!
! Compute the RMSD and the corresponding rotation matrix:
!
call s_rmsd_fast_calc_rmsd_and_rotation(A,-1,ncount,E0,rmsd,rot,ireturn)
!
! Nullify the coordinates of the responding geometric center vectors:
!
center1=0.0
center2=0.0
!
! Compute the accumulation required to obtain the coordinates of
! the geometric center vectors:
!
do i=1,num_atoms
   !
   if (do_count(i)) then
      !
      center1=center1+x1(:,i)
      center2=center2+x2(:,i)
      !
   end if
   !
end do
!
! Do average the accumulation to obtain the vectors of the centers:
!
center1=center1/ncount
center2=center2/ncount
!
! Use the rotation matrix computed in attempt to minimize the RMSD,
! to superimpose the coordinates of the compared molecules:
!
do i=1,num_atoms
   !
   ! Shift the atom coordinates of the first molecule with respect
   ! to their own geometric center vector:
   !
   x2(:,i)=x2(:,i)-center2
   !
   ! Backup the centered coordinates before starting the rotation.
   ! They will be subject of alteration during the rotation process
   ! that follows bellow!
   !
   x_tmp=x2(:,i)
   !
   ! Do the rotation direction by direction:
   !
   do j=1,3
      !
      x2(j,i)=rot(t(1,j))*x_tmp(t(2,j))+&
              rot(t(3,j))*x_tmp(t(4,j))+&
              rot(t(5,j))*x_tmp(t(6,j))
      !
   end do
   !
   ! Shift the coordinates of the atoms of the second molecule
   ! with respect to the geometric center vector of the atom coordinates
   ! of the first molecule:
   !
   x2(:,i)=x2(:,i)+center1
   !
end do

end subroutine s_rmsd_superimpose
