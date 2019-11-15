module postproc
  contains

  subroutine read_fortran_binary(fname,vartype,icount,string,iarray,rarray,ios)
  !------------------------------------------------------------------------
  !read unformatted fortran binary file
  !Input
  !    1)fname: file name
  !    2)vartype: (0:string; 1: integer; 2: real(8))
  !    3)icount:  number of values ( for string, icount must be 1)
  !
  !Output
  !    1)string: string outputs 
  !    2)iarray: integer array
  !    3)rarray: float array
  !    4)ios: (0: normal; <0: reaching the end) 
  !------------------------------------------------------------------------
    implicit none
    character(*),intent(in) :: fname
    integer,intent(in) :: vartype,icount
    integer,intent(out) :: ios
    character(len=20),intent(out) :: string 
    integer,intent(out) :: iarray(icount)
    real(kind=8),intent(out) :: rarray(icount)

    !local variable
    integer :: i
    logical, save :: lopen

    inquire(11,opened=lopen) 
    if(.not.lopen) open(11,file=fname,form='unformatted',status='old')
    
    if(vartype==0) read(11,iostat=ios)string
    if(vartype==1) read(11,iostat=ios)(iarray(i),i=1,icount)
    if(vartype==2) read(11,iostat=ios)(rarray(i),i=1,icount)
    if(ios>0) then
      stop 'wrong in outputs'
    elseif(ios<0) then
      close(11)
      return
    endif

  end subroutine read_fortran_binary

end module postproc
