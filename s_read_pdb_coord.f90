subroutine s_read_pdb_coord(pdb_file,num_atoms,coord)
!
! Reads the atomic coordinates declared inside the ATOM and HETATM
! records in PDB structure. Note, that one need to call first the
! subroutine "s_count_pdb_atoms" to estimate the number of those
! records.
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
character(len=*),intent(in) :: pdb_file
integer(C_INT),intent(in)   :: num_atoms
real(C_FLOAT),intent(out)   :: coord(3,num_atoms)
!
! Local variables:
!
integer(C_INT)              :: i
integer(C_INT)              :: stat
integer(C_INT)              :: io_stat
integer(C_INT),parameter    :: str_len=54
character(len=str_len)      :: buffer

open(666,file=pdb_file,iostat=stat)
!
i=0
!
do while(stat.eq.0)
   !
   read(666,'(A54)',iostat=stat) buffer
   !
   if (stat==0) then
      !
      if (buffer(1:6).eq.'ATOM  ') then
         !
         i=i+1
         !
         read(buffer(31:54),*,iostat=io_stat) coord(:,i)
         !
         if (io_stat/=0) then
            !
            print *
            print *,'CANNOT READ THE COORDINATES FROM THE INPUT PDB FILE!'
            print *
            !
            call exit(1)
            !
         end if
         !
      end if
      !
   end if
   !
end do
!
close(666)

end subroutine s_read_pdb_coord

