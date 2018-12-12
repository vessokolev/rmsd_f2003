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
real(C_FLOAT)              :: Sxx
real(C_FLOAT)              :: Sxy
real(C_FLOAT)              :: Sxz
real(C_FLOAT)              :: Syx
real(C_FLOAT)              :: Syy
real(C_FLOAT)              :: Syz
real(C_FLOAT)              :: Szx
real(C_FLOAT)              :: Szy
real(C_FLOAT)              :: Szz
real(C_FLOAT)              :: Szz2
real(C_FLOAT)              :: Syy2
real(C_FLOAT)              :: Sxx2
real(C_FLOAT)              :: Sxy2
real(C_FLOAT)              :: Syz2
real(C_FLOAT)              :: Sxz2
real(C_FLOAT)              :: Syx2
real(C_FLOAT)              :: Szy2
real(C_FLOAT)              :: Szx2
real(C_FLOAT)              :: SyzSzymSyySzz2
real(C_FLOAT)              :: Sxx2Syy2Szz2Syz2Szy2
real(C_FLOAT)              :: Sxy2Sxz2Syx2Szx2
real(C_FLOAT)              :: SxzpSzx
real(C_FLOAT)              :: SyzpSzy
real(C_FLOAT)              :: SxypSyx
real(C_FLOAT)              :: SyzmSzy
real(C_FLOAT)              :: SxzmSzx
real(C_FLOAT)              :: SxymSyx
real(C_FLOAT)              :: SxxpSyy
real(C_FLOAT)              :: SxxmSyy
real(C_FLOAT)              :: C(3)
real(C_FLOAT)              :: mxEigenV
real(C_FLOAT)              :: oldg
real(C_FLOAT)              :: b
real(C_FLOAT)              :: aa
real(C_FLOAT)              :: delta
real(C_FLOAT)              :: rms
real(C_FLOAT)              :: qsqr
real(C_FLOAT)              :: q1
real(C_FLOAT)              :: q2
real(C_FLOAT)              :: q3
real(C_FLOAT)              :: q4
real(C_FLOAT)              :: normq
real(C_FLOAT)              :: a11
real(C_FLOAT)              :: a12
real(C_FLOAT)              :: a13
real(C_FLOAT)              :: a14
real(C_FLOAT)              :: a21
real(C_FLOAT)              :: a22
real(C_FLOAT)              :: a23
real(C_FLOAT)              :: a24
real(C_FLOAT)              :: a31
real(C_FLOAT)              :: a32
real(C_FLOAT)              :: a33
real(C_FLOAT)              :: a34
real(C_FLOAT)              :: a41
real(C_FLOAT)              :: a42
real(C_FLOAT)              :: a43
real(C_FLOAT)              :: a44
real(C_FLOAT)              :: a2
real(C_FLOAT)              :: x2
real(C_FLOAT)              :: y2
real(C_FLOAT)              :: z2
real(C_FLOAT)              :: xy
real(C_FLOAT)              :: az
real(C_FLOAT)              :: zx
real(C_FLOAT)              :: ay
real(C_FLOAT)              :: yz
real(C_FLOAT)              :: ax
real(C_FLOAT)              :: a3344_4334
real(C_FLOAT)              :: a3244_4234
real(C_FLOAT)              :: a3243_4233
real(C_FLOAT)              :: a3143_4133
real(C_FLOAT)              :: a3144_4134
real(C_FLOAT)              :: a3142_4132
real(C_FLOAT),parameter    :: evecprec=1.0e-6
real(C_FLOAT),parameter    :: evalprec=1.0e-11
real(C_FLOAT)              :: a1324_1423
real(C_FLOAT)              :: a1224_1422
real(C_FLOAT)              :: a1223_1322
real(C_FLOAT)              :: a1124_1421
real(C_FLOAT)              :: a1123_1321
real(C_FLOAT)              :: a1122_1221

Sxx=A(1)
Sxy=A(2)
Sxz=A(3)
Syx=A(4)
Syy=A(5)
Syz=A(6)
Szx=A(7)
Szy=A(8)
Szz=A(9)
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
mxEigenV = E0
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
rms=sqrt(abs(2*(E0-mxEigenV)/npro))
rmsd=rms
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
rot(1)=a2+x2-y2-z2
rot(2)=2*(xy+az)
rot(3)=2*(zx-ay)
!
rot(4)=2*(xy-az)
rot(5)=a2-x2+y2-z2
rot(6)=2*(yz+ax)
!
rot(7)=2*(zx+ay)
rot(8)=2*(yz-ax)
rot(9)=a2-x2-y2+z2
!
ireturn=1

end subroutine s_rmsd_fast_calc_rmsd_and_rotation

