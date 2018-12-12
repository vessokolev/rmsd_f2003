subroutine s_rmsd_inner_product(num_atoms,do_count,x1,x2,A,E0)
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
! Version: 2018121200
!
use iso_c_binding,only:C_INT,C_FLOAT
!
implicit none
!
! Interface variables:
!
integer(C_INT),intent(in) :: num_atoms
logical,intent(in)        :: do_count(num_atoms)
real(C_FLOAT),intent(in)  :: x1(3,num_atoms)
real(C_FLOAT),intent(in)  :: x2(3,num_atoms)
real(C_FLOAT),intent(out) :: A(9)
real(C_FLOAT),intent(out) :: E0
!
! Local variables:
!
integer(C_INT)            :: i
integer(C_INT)            :: j
integer(C_INT)            :: k(2)
real(C_FLOAT)             :: G
real(C_FLOAT)             :: tmp

! Nullify the elements of the array 'A':
!
A=0.0
!
! ... and do the same for the value of 'G':
!
G=0.0
!
do i=1,num_atoms
   ! 
   if (do_count(i)) then
      !
      tmp=0.0
      !
      do j=1,3
         !
         k(1)=3*(j-1)
         !
         k(1)=k(1)+1
         k(2)=k(1)+2
         !
         A(k(1):k(2))=A(k(1):k(2))+x1(j,i)*x2(:,i)
         !
         tmp=tmp+x1(j,i)**2+x2(j,i)**2
         !
      end do
      !
      G=G+tmp
      !
   end if
   !
end do
!
E0=G/2

end subroutine s_rmsd_inner_product

