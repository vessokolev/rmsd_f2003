subroutine s_write_pdb_file(pdb_file,num_atoms,pdb_line,coord)
!
! Writes the atomic ATOM and HETATM records to a new PDB file.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2018121000
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT
!
implicit none
!
! Interface variables:
!
character(len=*),intent(in)  :: pdb_file
integer(C_INT),intent(in)    :: num_atoms
character(len=54),intent(in) :: pdb_line(num_atoms)
real(C_FLOAT),intent(in)     :: coord(3,num_atoms)
!
! Local variables:
!
integer(C_INT)               :: i
integer(C_INT)               :: stat
integer(C_INT),parameter     :: str_len=54

open(666,file=pdb_file,status='unknown',iostat=stat)
!
do i=1,num_atoms
   !
   write(666,fmt='(A30,3F8.3)') pdb_line(i)(1:30),coord(:,i)
   !
end do
!
close(666)

end subroutine s_write_pdb_file

