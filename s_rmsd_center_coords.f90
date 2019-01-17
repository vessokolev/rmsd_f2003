subroutine s_rmsd_center_coords(num_atoms,do_count,coord,ncount)
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
!    determination of the optimal rotational matrix for do_counted  superpositions."
!    Journal of Computational Chemistry 31(7):1561-1563
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! License: BSD
! Version: 2019011700
!
use iso_c_binding,only:C_INT,C_FLOAT
!
implicit none
!
! Interface variables:
!
integer(C_INT),intent(in)   :: num_atoms
logical,intent(in)          :: do_count(num_atoms)
real(C_FLOAT),intent(inout) :: coord(3,num_atoms)
integer(C_INT),intent(out)  :: ncount
!
! Local variables:
!
integer(C_INT)              :: i
real(C_FLOAT)               :: center(3)

! Nullify the components of the vector 'center':
!
center=0.0
!
! Initialize ncount
!
ncount=0
!
! Accumulate x, y, and z-coordinates of the atoms into the
! components of the vector 'center'. Count only the atoms for which
! 'do_count' component is set '.true.':
!
do i=1,num_atoms
   !
   if (do_count(i)) then
      !
      center=center+coord(:,i)
      !
      ncount=ncount+1
      !
   end if
   !
end do
!
! Average the accumulation and obtaint the coordinates of the
! geometric center vector:
!
center=center/ncount
!
! Shift the coordinates of the atoms with respect to the geometric center,
! which vector is estimated above:
!
do i=1,num_atoms
   !
   if (do_count(i)) then
      !
      coord(:,i)=coord(:,i)-center
      !
   end if
   !
end do

end subroutine s_rmsd_center_coords

