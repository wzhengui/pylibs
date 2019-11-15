module grid_info

  integer,parameter :: rkind=8,max_i34=6
  
  !Element geometry data
  integer,save :: ne
  integer,save,allocatable :: elnode(:,:)
  integer,save,allocatable :: elside(:,:)
  integer,save,allocatable :: ic3(:,:)
  integer,save,allocatable :: i34(:)
  real(rkind),save,allocatable :: xctr(:),yctr(:),dpe(:)
  real(rkind),save,allocatable :: area(:)

  !Node geometry data
  integer,save :: np
  integer,save :: mnei
  integer,save :: nx(max_i34-1,max_i34,max_i34)
  integer,save,allocatable :: nne(:),indel(:,:)
  real(rkind),save,allocatable :: xnd(:),ynd(:),dp(:)
  real(rkind),save,allocatable :: lon(:),lat(:),lonctr(:),latctr(:)

  !Side geometry data
  integer,save :: ns
  integer,save,allocatable :: isdel(:,:)
  integer,save,allocatable :: isidenode(:,:)
  real(rkind),save,allocatable :: xcj(:),ycj(:),dps(:)
  real(rkind),save,allocatable :: distj(:)

  !open and land boundary data
  integer,save :: nope,nland
  integer,save :: neta,nvel
  integer,save :: mnond,mnlnd,mnosd,mnoe
  integer,save,allocatable :: isbnd(:),isbnode(:,:),isbs(:,:),isbe(:,:) 
  !different definition from isbnd_global in SCHISM
  integer,save,allocatable :: nond(:),iond(:,:) !node in open bnd
  integer,save,allocatable :: nlnd(:),ilnd(:,:) !node in land bnd
  integer,save,allocatable :: nosd(:),iosd(:,:) !side in open bnd
  integer,save,allocatable :: noe(:),ioe(:,:)   !element in open bnd
 
end module grid_info

subroutine read_grid_info(Grid)
!--------------------------------------------------------------------------------
!this subrontine aquire the horizontal grid information
!construct element/nodes/sides relation from hgrid.gr3
!--------------------------------------------------------------------------------
  use grid_info
  implicit none
  character(*),intent(in):: Grid
  real(rkind) :: signa

  !----local variables------------------------------
  integer :: i,j,k,m,n,l,ip,ie,je,jsj
  integer :: itmp,itmp1,itmp2,n1,n2,nt,nn,n3,n4,jj,j2,j3,istat
  real(rkind) :: tmp,tmp1,tmp2,x1,x2,x3,x4,y1,y2,y3,y4
  character(len=200) :: errmsg
  logical :: found,lexist

  !calculate cyclic node index
  do k=3,max_i34 !type of poly
    do i=1,k
      do j=1,k-1
        nx(j,i,k)=i+j
        if(nx(j,i,k)>k) nx(j,i,k)=nx(j,i,k)-k
        if(nx(j,i,k)<1.or.nx(j,i,k)>k) then
          write(errmsg,*)'INIT: nx wrong',i,j,k,nx(j,i,k)
          call abort(errmsg)
        endif
      enddo !j
    enddo !i
  enddo !k


  !open grid file
  open(14,file=trim(adjustl(Grid)),status='old',iostat=istat)
  if(istat/=0) call abort('failed in open hgrid.gr3')
  
  read(14,*); read(14,*) ne,np

  !read lon lat information
  lexist=.false.; inquire(file='hgrid.ll',exist=lexist)
  if(lexist) then
    allocate(lon(np),lat(np),lonctr(ne),latctr(ne),stat=istat)
    if(istat/=0) call abort('failed in alloc. lon')
 
     open(15,file='hgrid.ll',status='old',iostat=istat)
     if(istat/=0) call abort('failed in open hgrid.ll')
     read(15,*); read(15,*)
     do i=1,np
       read(15,*) itmp,lon(i),lat(i),tmp
     enddo
     close(15)
     lonctr=0.0; latctr=0.0
   endif
  
  !read node coordinates, element to node table
  allocate(xnd(np),ynd(np),dp(np),xctr(ne),yctr(ne),dpe(ne),i34(ne),elnode(max_i34,ne),area(ne),stat=istat)
  if(istat/=0) call abort('failed in alloc. xnd')
  do i=1,np
    read(14,*) itmp,xnd(i),ynd(i),dp(i)
  enddo

  xctr=0.0; yctr=0.0; area=0.0
  do i=1,ne
    read(14,*) ie,i34(ie),(elnode(j,ie),j=1,i34(ie))
    do j=1,i34(ie)
       xctr(ie)=xctr(ie)+xnd(elnode(j,ie))/i34(ie)
       yctr(ie)=yctr(ie)+ynd(elnode(j,ie))/i34(ie)
       if(lexist) then
         lonctr(ie)=lonctr(ie)+lon(elnode(j,ie))/i34(ie)
         latctr(ie)=latctr(ie)+lat(elnode(j,ie))/i34(ie)
       endif
    enddo !j
    
    do j=1,i34(ie)-2
      j2=nx(1,j,i34(ie)); j3=nx(2,j,i34(ie))
      x1=xnd(elnode(1,ie));x2=xnd(elnode(j2,ie));x3=xnd(elnode(j3,ie));
      y1=ynd(elnode(1,ie));y2=ynd(elnode(j2,ie));y3=ynd(elnode(j3,ie));
      area(ie)=area(ie)+signa(x1,x2,x3,y1,y2,y3)
    enddo

    !dpe(i)=sum(dp(elnode(1:i34(i),i)))/i34(i)
    
    if(area(ie)<0) then
      write(errmsg,*)'area(ie)<0,wrong orientation:',ie,elnode(1:i34(ie),ie)
      call abort(errmsg)
    endif
  enddo !i


  !count number of element connected to each node
  allocate(nne(np),stat=istat)
  if(istat/=0) call abort('failed in alloc. nne')

  nne=0
  do ie=1,ne
    do i=1,i34(ie)
      ip=elnode(i,ie)
      nne(ip)=nne(ip)+1
    enddo !i
  enddo !ie

  !check hanging node
  found=.false.
  do ip=1,np
    if(nne(ip)==0) then
      write(11,*)'hanging node:',ip
      found=.true.
    endif
    if(found) call abort('hanging node found, check hgrid.gr3')
  enddo

  !maximum number of elements around a node
  mnei=0
  do ip=1,np
    mnei=max(mnei,nne(ip))
  enddo
  
  !buld node to-element table
  allocate(indel(mnei,np),stat=istat) 
  if(istat/=0) call abort('failed in alloc. indel')

  nne=0; indel=0
  do ie=1,ne
    do i=1,i34(ie)
      ip=elnode(i,ie)
      nne(ip)=nne(ip)+1
      indel(nne(ip),ip)=ie
    enddo !i
  enddo !ie
  
!check indel
!  do ip=1,np
!    do i=1,nne(ip)
!      write(16,*)'node-elem',ip,nne(ip),i,indel(i,ip)
!    enddo
!  enddo

  !build element-side-element table
  allocate(ic3(max_i34,ne),stat=istat) 
  if(istat/=0) call abort('failed in alloc. indel')
  
  do ie=1,ne
    do i=1,i34(ie)
      ic3(i,ie)=0 !index for boundary sides
      n1=elnode(nx(1,i,i34(ie)),ie)
      n2=elnode(nx(2,i,i34(ie)),ie)
      jj=0 !count
      do k=1,nne(n1)  !search around node n1
        je=indel(k,n1)
        if(je/=ie) then
           do m=1,i34(je)
             if(elnode(m,je)==n2) then
               jj=jj+1
               ic3(i,ie)=je; exit
             endif
           enddo !m
        endif
      enddo! k

      if(jj>1) then
         write(errmsg,*)'more than one neigbr issue,',ie,n1,n2,jj
         call abort(errmsg)
      endif

      je=ic3(i,ie)
      if(je/=0) then !check whether nodes have opposite orientation
        do k=1,i34(je)
          if(elnode(nx(1,k,i34(je)),je)==n1.and.elnode(nx(2,k,i34(je)),je)==n2) then
            write(*,*)'element ',ie, ' and ',je, ' have opposite orientation',n1,n2
            call abort('Elements have opposite orientation ')
          endif
        enddo
      endif

    enddo !i
  enddo !ie
  
  !build element-side index table
  allocate(elside(max_i34,ne),stat=istat)
  if(istat/=0) call abort('failed in alloc. elside')

  ns=0;elside=0
  do ie=1,ne
    do j=1,i34(ie)
      if(ic3(j,ie)==0.or.ie<ic3(j,ie)) then !found a new side
        ns=ns+1
        elside(j,ie)=ns
        if(ic3(j,ie)/=0) then !add side index for neighbor element
          je=ic3(j,ie)
          l=0
          do k=1,i34(je)
            if(ic3(k,je)==ie) then
              l=k
              exit
            endif
          enddo
          if(l==0) then
            write(errmsg,'(a,10i6)')'grid_subs: wrong ball info',ie,j,ns 
            call abort(errmsg)
          endif
          elside(l,je)=ns
        endif
      endif !ic3(j,ie)==0.and.ie<ic3(j,ie)
    enddo
  enddo !ie

  !check elside
  found=.false.
  do ie=1,ne
    do i=1,i34(ie)
      if(elside(i,ie)==0) then
        write(11,'(a,2i6)')'side is missing ',ie,i
        found=.true.
      endif
    enddo
  enddo
  if(found) call abort('side is missing') 
  
  if(ns<ne.or.ns<np) then !side number should be larger than elem. and node num.
    call abort('ns<ne or ns<np')
  endif

  !alloc. side arrays
  allocate(isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns),dps(ns),distj(ns),stat=istat)
  if(istat/=0) call abort('failed in alloc. isdel')


  !side-element, side-node table
  dpe=-huge(1.0)
  do ie=1,ne
    do j=1,i34(ie)
      jsj=elside(j,ie)
      n1=elnode(nx(1,j,i34(ie)),ie)
      n2=elnode(nx(2,j,i34(ie)),ie)
      if(ic3(j,ie)==0.or.ie<ic3(j,ie)) then !new side found
        isdel(1,jsj)=ie
        isdel(2,jsj)=ic3(j,ie)
        isidenode(1,jsj)=n1
        isidenode(2,jsj)=n2
      else  !ie>ic3(j,ie)
        isdel(1,jsj)=ic3(j,ie)
        isdel(2,jsj)=ie
        isidenode(1,jsj)=n2
        isidenode(2,jsj)=n1
      endif
      xcj(jsj)=(xnd(n1)+xnd(n2))/2.0d0
      ycj(jsj)=(ynd(n1)+ynd(n2))/2.0d0
      dps(jsj)=(dp(n1)+dp(n2))/2.0d0
      distj(jsj)=sqrt((xnd(n2)-xnd(n1))**2+(ynd(n2)-ynd(n1))**2)
      if(dps(jsj)>dpe(ie)) dpe(ie)=dps(jsj)
    enddo
  enddo !ie

!--------------------this part read bnd.------------------------------------
  !read open boundary info.
  read(14,*) nope; read(14,*) neta

  mnond=0 !max number of nodes per segment
  nt=0   !total open bnd. node
  do k=1,nope
    read(14,*) nn
    mnond=max(mnond,nn)
    nt=nt+nn
    do i=1,nn; read(14,*); enddo
  enddo
  if(neta/=nt) then
    write(errmsg,*) 'neta/= total # of open bnd nodes',neta,nt
    call abort(errmsg)
  endif
  
  !alloc. arrays
  allocate(isbnd(np),nond(nope),iond(mnond,nope),stat=istat)
  if(istat/=0) call abort('failed in alloc. isbnd')
  !read open bnd. segments and nodes
  rewind(14); read(14,*); read(14,*)
  do i=1,np; read(14,*); enddo
  do i=1,ne; read(14,*); enddo
  read(14,*); read(14,*)
  isbnd=0; nond=0; iond=0
  do k=1,nope
    read(14,*) nn
    do i=1,nn
      read(14,*)ip
      nond(k)=nond(k)+1
      iond(nond(k),k)=ip
      isbnd(ip)=1
    enddo !i
    if(iond(1,k)==iond(nond(k),k)) then
      write(errmsg,*) 'Looped open bnd:',k
      call abort(errmsg)
    endif
  enddo !k
  
  !aquire land bnd. seg and nodes
  rewind(14); read(14,*); read(14,*)
  do i=1,np; read(14,*); enddo
  do i=1,ne; read(14,*); enddo
  read(14,*); read(14,*)
  do k=1,nope; read(14,*) nn; do i=1,nn; read(14,*); enddo; enddo
  read(14,*) nland
  read(14,*) nvel

  mnlnd=0 !max number of nodes per segment
  nt=0   !total land bnd. node
  do k=1,nland
    read(14,*) nn
    mnlnd=max(mnlnd,nn)
    nt=nt+nn
    do i=1,nn; read(14,*); enddo
  enddo !k
  if(nvel/=nt) then
    write(errmsg,*) 'nvel/= total # of land bnd nodes',nvel,nt
    call abort(errmsg)
  endif

  !alloc. arrayrs for land bnd. 
  allocate(nlnd(nland),ilnd(mnlnd,nland),stat=istat)
  if(istat/=0) call abort('failed in alloc. nlnd')
  rewind(14); read(14,*); read(14,*)
  do i=1,np; read(14,*); enddo
  do i=1,ne; read(14,*); enddo
  read(14,*); read(14,*)
  do k=1,nope; read(14,*) nn; do i=1,nn; read(14,*); enddo; enddo
  read(14,*); read(14,*)
  nlnd=0; ilnd=0
  do k=1,nland
    read(14,*)nn
    do i=1,nn
      read(14,*)ip
      nlnd(k)=nlnd(k)+1
      ilnd(nlnd(k),k)=ip
      isbnd(ip)=isbnd(ip)+2
    enddo !i
  enddo !k
  close(14)

  !open bnd node flag isbnode (different from isbnd)
  allocate(isbnode(2,np),stat=istat)
  if(istat/=0) call abort('failed in alloc. isbnode')
  isbnode=0
  do k=1,nope
    do j=1,nond(k)
      ip=iond(j,k)
      if(isbnode(1,ip)==0) then
        isbnode(1,ip)=k
      elseif(isbnode(2,ip)==0) then
        isbnode(2,ip)=k
      else
        call abort('wrong for isbnode')
      endif
    enddo
  enddo !k
 
  !Allocate and count global boundary side assigned to open boundary segments 
  if(allocated(nosd)) deallocate(nosd);
  allocate(nosd(nope),stat=istat)
  if(istat/=0) call abort('failed in alloc. nosd')
  nosd=0
  do ie=1,ne
    Lsg1: do j=1,i34(ie)
      if(ic3(j,ie)==0) then
         n1=elnode(nx(1,j,i34(ie)),ie)
         n2=elnode(nx(2,j,i34(ie)),ie)
         do k=1,nope
           if((isbnode(1,n1)==k.or.isbnode(2,n1)==k).and. &
             & (isbnode(1,n2)==k.or.isbnode(2,n2)==k)) then
             nosd(k)=nosd(k)+1
             cycle Lsg1
           endif 
         enddo !k
      endif !ic3
    enddo Lsg1 ! j
  enddo !ie

  mnosd=0
  do k=1,nope
    mnosd=max(mnosd,nosd(k))
  enddo
 
  !Allocate and count global boundary side assigned to open boundary segments 
  if(allocated(iosd)) deallocate(iosd);
  allocate(iosd(mnosd,nope),stat=istat)
  if(istat/=0) call abort('failed in alloc. iosd')
  iosd=0; nosd=0
  do ie=1,ne
    Lsg2: do j=1,i34(ie)
      if(ic3(j,ie)==0) then
         n1=elnode(nx(1,j,i34(ie)),ie)
         n2=elnode(nx(2,j,i34(ie)),ie)
         do k=1,nope
           if((isbnode(1,n1)==k.or.isbnode(2,n1)==k).and. &
             & (isbnode(1,n2)==k.or.isbnode(2,n2)==k)) then
             nosd(k)=nosd(k)+1
             iosd(nosd(k),k)=elside(j,ie)
             cycle Lsg2
           endif 
         enddo !k
      endif !ic3
    enddo Lsg2 !j
  enddo !ie
  
  !Allocate and count boundary elements assigned to open boundary segments 
  if(allocated(noe)) deallocate(noe)
  allocate(noe(nope),stat=istat)
  if(istat/=0) call abort('failed in alloc. noe')
  noe=0
  Leg1: do ie=1,ne
    do j=1,i34(ie)
      ip=elnode(j,ie)
      do k=1,nope
        if(isbnode(1,ip)==k.or.isbnode(2,ip)==k) then
          noe(k)=noe(k)+1
          cycle Leg1
        endif
      enddo !k
    enddo !j
  enddo Leg1
  
  mnoe=0
  do k=1,nope
    mnoe=max(mnoe,noe(k))
  enddo
  
  !Allocate and count boundary elements assigned to open boundary segments 
  if(allocated(ioe)) deallocate(ioe)
  allocate(ioe(mnoe,nope),stat=istat)
  if(istat/=0) call abort('failed in alloc. noe')
  noe=0; ioe=0;
  Leg2: do ie=1,ne
    do j=1,i34(ie)
      ip=elnode(j,ie)
      do k=1,nope
        if(isbnode(1,ip)==k.or.isbnode(2,ip)==k) then
          noe(k)=noe(k)+1
          ioe(noe(k),k)=ie
          cycle Leg2
        endif
      enddo 
    enddo
  enddo Leg2

  !Allocate and classify global boundary elements and sides
  allocate(isbe(2,ne),isbs(2,ns),stat=istat)
  if(istat/=0) call abort('failed in alloc. isbe')
  isbe=0; isbs=0
  do k=1,nope
    do j=1,noe(k)
      isbe(1,ioe(j,k))=k
      isbe(2,ioe(j,k))=j
    enddo
    do j=1,nosd(k)
      isbs(1,iosd(j,k))=k 
      isbs(2,iosd(j,k))=j
    enddo
  enddo !k
  do ie=1,ne
    do j=1,i34(ie)
      jsj=elnode(j,ie) 
      if(ic3(j,ie)==0.and.isbs(1,jsj)==0) isbs(1,jsj)=-1
    enddo !j
  enddo !ie
   
 
end subroutine read_grid_info

subroutine abort(errmsg)
  implicit none
  character(*),optional,intent(in) :: errmsg

  if(present(errmsg)) then
    write(*,*)errmsg
  endif
  stop
end subroutine abort

function signa(x1,x2,x3,y1,y2,y3)
!-------------------------------------------------------------------------------
! Compute signed area formed by pts 1,2,3
!-------------------------------------------------------------------------------
  use grid_info, only : rkind
  implicit none
  real(rkind) :: signa
  real(rkind),intent(in) :: x1,x2,x3,y1,y2,y3

  signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2d0

end function signa

subroutine write_obe
  use grid_info
  implicit none

  !local variables
  integer :: i,j
  real(kind=rkind) :: xtmp, ytmp
  logical :: lexist

!---------------------------------------------------------------
!write obe.out
!---------------------------------------------------------------
  !element 
  open(32,file='obe.out')
  write(32,*) nope,' # of open bnd'
  write(32,*) 'Element list:'
  do i=1,nope
    write(32,*) noe(i),' bnd #',i
    do j=1,noe(i)
       write(32,*) j,ioe(j,i)
    enddo !j
  enddo !i

  !side
  write(32,*) 'Side list:'
  do i=1,nope
    write(32,*) nosd(i),' bnd #',i
    do j=1,nosd(i)
        write(32,*) j,iosd(j,i)
    enddo
  enddo !i
  close(32)

!---------------------------------------------------------------
!sidecenters.bp
!---------------------------------------------------------------
  open(32,file='sidecenters.bp')
  write(32,*) 'Sidegrid'
  write(32,*) ns
  do i=1,ns
      write(32,*) i,xcj(i),ycj(i),real(dps(i))
  enddo !i
  close(32)

!---------------------------------------------------------------
!centers.bp
!---------------------------------------------------------------
  open(32,file='centers.bp',status='replace')
  write(32,*) 'centers pts'
  write(32,*) ne
  do i=1,ne
    write(32,*) i,xctr(i),yctr(i),real(dpe(i))
  enddo !i
  close(32)

!---------------------------------------------------------------
!centers.ll if hgrid.ll exists
!---------------------------------------------------------------
  inquire(file='hgrid.ll',exist=lexist)
  if(lexist) then
    open(32,file='centers.ll',status='replace')
    write(32,*) 'centers pts in lon/lat'
    write(32,*) ne
    do i=1,ne
      xtmp=sum(lon(elnode(1:i34(i),i)))/i34(i)
      ytmp=sum(lat(elnode(1:i34(i),i)))/i34(i)
      write(32,*) i,xtmp,ytmp,real(dpe(i))
    enddo !i
    close(32)
  endif !lexist
end subroutine write_obe


