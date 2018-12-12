program demo
!
! This program illustrates how does the QCP method for computing RMSD
! between two molecules work. The first molecule (stored in 'pdb_1') is
! set a reference. Then second molecule (stored in 'pdb_2') will be
! translated and rotated until matches, as much as possible, the referent
! molecule. The result of the translation and rotation will be stored
! in 'pdb_3' file. One can preview the inputs and output by using VMD
! or UCSF Chimera.
!
! To cite QCP:
!
! 1. Douglas L. Theobald (2005) "Rapid calculation of RMSD using a
! quaternion-based
!    characteristic polynomial." Acta Crystallographica A 61(4):478-480.
!
! 2. Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2010) "Fast
!    determination of the optimal rotational matrix for do_counted
!    superpositions."
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
integer(C_INT)                :: num_atoms_1
integer(C_INT)                :: num_atoms_2
real(C_FLOAT),allocatable     :: x1(:,:)
real(C_FLOAT),allocatable     :: x2(:,:)
real(C_FLOAT)                 :: rmsd
logical,allocatable           :: do_count(:)
character(len=54),allocatable :: pdb_line(:)
character(len=4096)           :: pdb_1
character(len=4096)           :: pdb_2
character(len=4096)           :: pdb_3

! The referent structure:
!
pdb_1='molecule_1.pdb'
!
! The molecule to compare with 'molecule_1':
!
pdb_2='molecule_2.pdb'
!
! The rotated version of the 'molecule 2':
!
pdb_3='molecule_2_rotated.pdb'
!
call s_count_pdb_atoms(pdb_1,num_atoms_1)
!
allocate(x1(3,num_atoms_1))
!
allocate(pdb_line(num_atoms_1))
!
call s_read_pdb_coord(pdb_1,num_atoms_1,x1)
!
call s_read_pdb_atom_lines(pdb_1,num_atoms_1,pdb_line)
!
call s_count_pdb_atoms(pdb_2,num_atoms_2)
!
if (num_atoms_1/=num_atoms_2) then
   !
   print *
   print *,'This program does not support the case when'
   print *,'the compared molecules have different number'
   print *,'of atoms.'
   print *
   !
   call exit(1)
   !
end if
!
allocate(x2(3,num_atoms_2))
!
call s_read_pdb_coord(pdb_2,num_atoms_2,x2)
!
allocate(do_count(ubound(x1,2)))
!
! Compare the molecules by involving all atoms:
!
do_count=.true.

call s_rmsd_superimpose(num_atoms_1,do_count,x1,x2,rmsd)
!
print *
print *,'The computed RMSD between the two molecules is: ',&
        rmsd,'Angstroms.'
print *
print *,'The coordinates of the second molecule, '//&
        'superimposed to those of the first one, are '//&
        'stored inside the file ',trim(adjustl(pdb_3))
print *
!
call s_write_pdb_file(pdb_3,num_atoms_1,pdb_line,x2)

end program demo
