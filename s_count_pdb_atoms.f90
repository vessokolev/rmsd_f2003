subroutine s_count_pdb_atoms(pdb_file,num_atoms)
!
! Reads the number of ATOM or HETATM records declared in PDB file.
! Single model hosted inside the PDB file is assumed.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2018121000
! License: GPLv2
!
use iso_c_binding,only:C_INT
!
implicit none
!
! Interface variables:
!
character(len=*),intent(in) :: pdb_file
integer(C_INT),intent(out)  :: num_atoms
!
! Local variables:
!
integer(C_INT)              :: stat
integer(C_INT),parameter    :: str_len=54
character(len=str_len)      :: buffer

open(666,file=pdb_file,iostat=stat)
!
num_atoms=0
!
do while(stat.eq.0)
   !
   read(666,'(A54)',iostat=stat) buffer
   !
   if (stat==0) then
      !
      if (buffer(1:6).eq.'ATOM  ') then
         !
         num_atoms=num_atoms+1
         !
      end if
      !
   end if
   !
end do
!
close(666)

end subroutine s_count_pdb_atoms

