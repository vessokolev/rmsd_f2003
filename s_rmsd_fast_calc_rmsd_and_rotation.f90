subroutine s_rmsd_fast_calc_rmsd_and_rotation(A,minScore,npro,E0,rmsd,rot,ireturn)
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
! NOTE: Use double precision in the computations presented bellow. Eploying
!       single precision floating point data can turn the elements of the
!       rotation matrix "rot" into zeroes. That will result into total
!       failure of the coordinate super imposition.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! License: BSD
! Version: 2019011500
!
use iso_c_binding,only:C_INT,C_FLOAT,C_DOUBLE
!
implicit none 
!
! Interface variables:
!
real(C_FLOAT),intent(in)   ::  A(9)
integer(C_INT),intent(in)  :: minScore
integer(C_INT),intent(in)  :: npro
real(C_FLOAT),intent(in)   :: E0
real(C_FLOAT),intent(out)  :: rmsd
real(C_FLOAT),intent(out)  :: rot(9)
integer(C_INT),intent(out) :: ireturn
!
! Local variables:
!
integer(C_INT)             :: i
real(C_DOUBLE)             :: Sxx
real(C_DOUBLE)             :: Sxy
real(C_DOUBLE)             :: Sxz
real(C_DOUBLE)             :: Syx
real(C_DOUBLE)             :: Syy
real(C_DOUBLE)             :: Syz
real(C_DOUBLE)             :: Szx
real(C_DOUBLE)             :: Szy
real(C_DOUBLE)             :: Szz
real(C_DOUBLE)             :: Szz2
real(C_DOUBLE)             :: Syy2
real(C_DOUBLE)             :: Sxx2
real(C_DOUBLE)             :: Sxy2
real(C_DOUBLE)             :: Syz2
real(C_DOUBLE)             :: Sxz2
real(C_DOUBLE)             :: Syx2
real(C_DOUBLE)             :: Szy2
real(C_DOUBLE)             :: Szx2
real(C_DOUBLE)             :: SyzSzymSyySzz2
real(C_DOUBLE)             :: Sxx2Syy2Szz2Syz2Szy2
real(C_DOUBLE)             :: Sxy2Sxz2Syx2Szx2
real(C_DOUBLE)             :: SxzpSzx
real(C_DOUBLE)             :: SyzpSzy
real(C_DOUBLE)             :: SxypSyx
real(C_DOUBLE)             :: SyzmSzy
real(C_DOUBLE)             :: SxzmSzx
real(C_DOUBLE)             :: SxymSyx
real(C_DOUBLE)             :: SxxpSyy
real(C_DOUBLE)             :: SxxmSyy
real(C_DOUBLE)             :: C(3)
real(C_DOUBLE)             :: mxEigenV
real(C_DOUBLE)             :: oldg
real(C_DOUBLE)             :: b
real(C_DOUBLE)             :: aa
real(C_DOUBLE)             :: delta
real(C_DOUBLE)             :: rms
real(C_DOUBLE)             :: qsqr
real(C_DOUBLE)             :: q1
real(C_DOUBLE)             :: q2
real(C_DOUBLE)             :: q3
real(C_DOUBLE)             :: q4
real(C_DOUBLE)             :: normq
real(C_DOUBLE)             :: a11
real(C_DOUBLE)             :: a12
real(C_DOUBLE)             :: a13
real(C_DOUBLE)             :: a14
real(C_DOUBLE)             :: a21
real(C_DOUBLE)             :: a22
real(C_DOUBLE)             :: a23
real(C_DOUBLE)             :: a24
real(C_DOUBLE)             :: a31
real(C_DOUBLE)             :: a32
real(C_DOUBLE)             :: a33
real(C_DOUBLE)             :: a34
real(C_DOUBLE)             :: a41
real(C_DOUBLE)             :: a42
real(C_DOUBLE)             :: a43
real(C_DOUBLE)             :: a44
real(C_DOUBLE)             :: a2
real(C_DOUBLE)             :: x2
real(C_DOUBLE)             :: y2
real(C_DOUBLE)             :: z2
real(C_DOUBLE)             :: xy
real(C_DOUBLE)             :: az
real(C_DOUBLE)             :: zx
real(C_DOUBLE)             :: ay
real(C_DOUBLE)             :: yz
real(C_DOUBLE)             :: ax
real(C_DOUBLE)             :: a1324_1423
real(C_DOUBLE)             :: a1224_1422
real(C_DOUBLE)             :: a1223_1322
real(C_DOUBLE)             :: a1124_1421
real(C_DOUBLE)             :: a1123_1321
real(C_DOUBLE)             :: a1122_1221
real(C_DOUBLE)             :: a3344_4334
real(C_DOUBLE)             :: a3244_4234
real(C_DOUBLE)             :: a3243_4233
real(C_DOUBLE)             :: a3143_4133
real(C_DOUBLE)             :: a3144_4134
real(C_DOUBLE)             :: a3142_4132
real(C_DOUBLE),parameter   :: evecprec=1.0e-6
real(C_DOUBLE),parameter   :: evalprec=1.0e-11

Sxx=dble(A(1))
Sxy=dble(A(2))
Sxz=dble(A(3))
Syx=dble(A(4))
Syy=dble(A(5))
Syz=dble(A(6))
Szx=dble(A(7))
Szy=dble(A(8))
Szz=dble(A(9))
!
Sxx2=Sxx*Sxx
Syy2=Syy*Syy
Szz2=Szz*Szz
!
Sxy2=Sxy*Sxy
Syz2=Syz*Syz
Sxz2=Sxz*Sxz
!
Syx2=Syx*Syx
Szy2=Szy*Szy
Szx2=Szx*Szx
!
SyzSzymSyySzz2=2*(Syz*Szy-Syy*Szz)
Sxx2Syy2Szz2Syz2Szy2=Syy2+Szz2-Sxx2+Syz2+Szy2
!
C(3)=-2*(Sxx2+Syy2+Szz2+Sxy2+Syx2+Sxz2+Szx2+Syz2+Szy2)
C(2)= 8*(Sxx*Syz*Szy+Syy*Szx*Sxz+Szz*Sxy*Syx-Sxx*Syy*Szz-&
         Syz*Szx*Sxy-Szy*Syx*Sxz)
!
SxzpSzx=Sxz+Szx
SyzpSzy=Syz+Szy
SxypSyx=Sxy+Syx
SyzmSzy=Syz-Szy
SxzmSzx=Sxz-Szx
SxymSyx=Sxy-Syx
SxxpSyy=Sxx+Syy
SxxmSyy=Sxx-Syy
!
Sxy2Sxz2Syx2Szx2=Sxy2+Sxz2-Syx2-Szx2
!
C(1)=Sxy2Sxz2Syx2Szx2*Sxy2Sxz2Syx2Szx2&
     +(Sxx2Syy2Szz2Syz2Szy2+SyzSzymSyySzz2)&
     *(Sxx2Syy2Szz2Syz2Szy2-SyzSzymSyySzz2)&
     +(-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz))&
     *(-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))&
     +(-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz))&
     *(-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))&
     +(+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz))&
     *(-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))&
     +(+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz))&
     *(-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz))
!
! Newton-Raphson
!
mxEigenV = dble(E0)
!
do i=1,50
   !
   oldg=mxEigenV
   x2=mxEigenV*mxEigenV
   b=(x2+C(3))*mxEigenV
   aa=b+C(2)
   delta=((aa*mxEigenV+C(1))/(2.0*x2*mxEigenV+b+aa))
   mxEigenV=mxEigenV-delta
   !
   if (abs(mxEigenV-oldg).lt.abs(evalprec*mxEigenV)) then
      exit
   end if
   !
end do
!
if (i.eq.50) then
   write(*,*) "More than",i,"iterations needed!"
endif
!
! the abs() is to guard against extremely small, but *negative* numbers due to floating point error
!
rms=sqrt(abs(2*(dble(E0)-mxEigenV)/npro))
rmsd=real(rms)
!
if ((minScore.gt.0).and.(rms.lt.minScore)) then
   ireturn=-1 ! Don't bother with rotation. 
   return
endif
!
a11=SxxpSyy+Szz-mxEigenV
a12=SyzmSzy
a13=-SxzmSzx
a14=SxymSyx
!
a21=SyzmSzy
a22=SxxmSyy-Szz-mxEigenV
a23=SxypSyx
a24=SxzpSzx
!
a31=a13
a32=a23
a33=Syy-Sxx-Szz-mxEigenV
a34=SyzpSzy
!
a41=a14
a42=a24
a43=a34
a44=Szz-SxxpSyy-mxEigenV
!
a3344_4334=a33*a44-a43*a34
a3244_4234=a32*a44-a42*a34
a3243_4233=a32*a43-a42*a33
a3143_4133=a31*a43-a41*a33
a3144_4134=a31*a44-a41*a34
a3142_4132=a31*a42-a41*a32
!
q1= a22*a3344_4334-a23*a3244_4234+a24*a3243_4233
q2=-a21*a3344_4334+a23*a3144_4134-a24*a3143_4133
q3= a21*a3244_4234-a22*a3144_4134+a24*a3142_4132
q4=-a21*a3243_4233+a22*a3143_4133-a23*a3142_4132
!
qsqr=q1*q1+q2*q2+q3*q3+q4*q4
!
! The following code tries to calculate another column in the adjoint matrix when the norm of the 
! current column is too small. Usually this block never becomes activated.
!
if (qsqr.lt.evecprec) then
   !
   q1= a12*a3344_4334-a13*a3244_4234+a14*a3243_4233
   q2=-a11*a3344_4334+a13*a3144_4134-a14*a3143_4133
   q3= a11*a3244_4234-a12*a3144_4134+a14*a3142_4132
   q4=-a11*a3243_4233+a12*a3143_4133-a13*a3142_4132
   !
   qsqr=q1*q1+q2*q2+q3*q3+q4*q4
   !
   if (qsqr.lt.evecprec) then
      !
      a1324_1423=a13*a24-a14*a23
      a1224_1422=a12*a24-a14*a22
      a1223_1322=a12*a23-a13*a22
      a1124_1421=a11*a24-a14*a21
      a1123_1321=a11*a23-a13*a21
      a1122_1221=a11*a22-a12*a21
      !
      q1= a42*a1324_1423-a43*a1224_1422+a44*a1223_1322
      q2=-a41*a1324_1423+a43*a1124_1421-a44*a1123_1321
      q3= a41*a1224_1422-a42*a1124_1421+a44*a1122_1221
      q4=-a41*a1223_1322+a42*a1123_1321-a43*a1122_1221
      !
      qsqr=q1*q1+q2*q2+q3*q3+q4*q4
      !
      if (qsqr.lt.evecprec) then
         !
         q1= a32*a1324_1423-a33*a1224_1422+a34*a1223_1322
         q2=-a31*a1324_1423+a33*a1124_1421-a34*a1123_1321
         q3= a31*a1224_1422-a32*a1124_1421+a34*a1122_1221
         q4=-a31*a1223_1322+a32*a1123_1321-a33*a1122_1221
         !
         qsqr=q1*q1+q2*q2+q3*q3+q4*q4
         !
         if (qsqr.lt.evecprec) then
            !
            ! if qsqr is still too small, return the identity matrix
            ! (that means no rotation of the atoms is reqiured)
            !
            rot(1)=1.0
            rot(5)=1.0
            rot(9)=1.0
            !
            rot(2)=0.0
            rot(3)=0.0
            rot(4)=0.0
            rot(6)=0.0
            rot(7)=0.0
            rot(8)=0.0
            !
            ireturn=0
            !
            return
            !
         end if
         !
      end if
      !
   end if
   !
end if
!
normq=sqrt(qsqr)
!
q1=q1/normq
q2=q2/normq
q3=q3/normq
q4=q4/normq
!
a2=q1*q1
x2=q2*q2
y2=q3*q3
z2=q4*q4
!
xy=q2*q3
az=q1*q4
zx=q4*q2
ay=q1*q3
yz=q3*q4
ax=q1*q2
!
rot(1)=real(a2+x2-y2-z2)
rot(2)=real(2*(xy+az))
rot(3)=real(2*(zx-ay))
!
rot(4)=real(2*(xy-az))
rot(5)=real(a2-x2+y2-z2)
rot(6)=real(2*(yz+ax))
!
rot(7)=real(2*(zx+ay))
rot(8)=real(2*(yz-ax))
rot(9)=real(a2-x2-y2+z2)
!
ireturn=1

end subroutine s_rmsd_fast_calc_rmsd_and_rotation

