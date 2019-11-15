!--------------------------------------------------------------
!written by ZG on 5/26/2016
!generate hotstart based on previous hotstart.in
!input: 
!     old grid: bg.gr3, vgrid.bg, hotstart.in.0 
!     new grid: fg.gr3, vgrid.fg
!output:
!     hotstart.in
!--------------------------------------------------------------
!pgf90 -O2 -mcmodel=medium  -Bstatic -o Interp_hotstart compute_zcor.f90 Interp_hotstart.F90

module var_glbl
  integer, parameter :: nbyte=4,rkind=8,rkind2=4,ntracers=2
  real(kind=rkind), parameter :: h0=0.1 
  integer :: nxq(3,4,4)

  !bg grid
  integer :: ne,np,ns,nvrt,ivcor
  integer, save,allocatable :: i34(:),elnode(:,:),nne(:),indel(:,:),elside(:,:), &
                             & ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),kbe(:)
  real(kind=rkind),save,allocatable :: xnd(:),ynd(:),dp(:),zcor(:,:),xctr(:),yctr(:)
  real(kind=rkind2),save,allocatable :: ztot(:),sigma(:),sigma_lcl(:,:)
  real(kind=rkind2),save :: kz,h_s,h_c,theta_b,theta_f 

  !fg grid
  integer :: ne_fg,np_fg,ns_fg,nvrt_fg,ivcor_fg
  integer, save,allocatable :: i34_fg(:),elnode_fg(:,:),nne_fg(:),indel_fg(:,:),elside_fg(:,:), &
                             & ic3_fg(:,:),isdel_fg(:,:),isidenode_fg(:,:),kbp_fg(:),kbe_fg(:)
  real(kind=rkind),save,allocatable :: xnd_fg(:),ynd_fg(:),dp_fg(:),zcor_fg(:,:),xctr_fg(:),yctr_fg(:)
  real(kind=rkind2),save,allocatable :: ztot_fg(:),sigma_fg(:),sigma_lcl_fg(:,:)
  real(kind=rkind2),save,allocatable :: kz_fg,h_s_fg,h_c_fg,theta_b_fg,theta_f_fg 

  !interpolation
  integer, save, allocatable :: ind2Dp(:),ind3Dp(:,:),ind2De(:),ind3De(:,:)

  !variables from hotstart.in.0
  real(kind=rkind),save,allocatable :: tr_el(:,:,:),tr_nd(:,:,:),tr_sed(:,:) 


  
end module var_glbl


program Interp_hotstart
!----------------------------------------------------------------
!interpolation 
!----------------------------------------------------------------
  use var_glbl

  !local variables
  implicit none
  integer :: i,j,k,l,m,n,ie,isum,istat,itmp,itmp1,itmp2,itmp3
  real(kind=rkind) :: dist,mdist,xi,yi,zi,xi_fg,yi_fg,zi_fg
  real(kind=rkind) :: rtmp,rtmp1,rtmp2,rtmp3,tmp

  !read bg grid info
  call read_grid
  !read fg grid info
  call read_grid_fg

  !alloc. variables  
  allocate(ind2Dp(np_fg),ind3Dp(nvrt_fg,np_fg),ind2De(ne_fg),ind3De(nvrt_fg,ne_fg), &
          & tr_el(nvrt,ntracers+1,ne),tr_nd(nvrt,2*ntracers+7,np),tr_sed(40,ne), stat=istat)
  if(istat/=0) stop 'failed in alloc ind2D'

  !interpolation index for node
  do i=1,np_fg
    l=0; mdist=1.d6;
    do j=1,np
      dist=sqrt((xnd_fg(i)-xnd(j))**2+(ynd_fg(i)-ynd(j))**2) 
      if(dist<mdist) then
        mdist=dist
        l=j
      endif !dist
    enddo !j
    if(l==0) stop 'can not find mdist'
    ind2Dp(i)=l

    do k=kbp_fg(i),nvrt_fg
      l=0;  mdist=1.d6
      do m=kbp(ind2Dp(i)),nvrt
        dist=abs(zcor_fg(k,i)-zcor(m,ind2Dp(i)))
        if(dist<mdist) then
          mdist=dist
          l=m
        endif !dist
      enddo  !m
      if(l==0) stop 'can not find mdist in zcor'
      ind3Dp(k,i)=l
      if(k==kbp_fg(i)) ind3Dp(1:kbp_fg(i)-1,i)=l
    enddo !k
  enddo !i

  !interpolation index for elem.
  do i=1,ne_fg
    l=0; mdist=1.d6;
    xi_fg=sum(xnd_fg(elnode_fg(1:i34_fg(i),i)))/i34_fg(i)
    yi_fg=sum(ynd_fg(elnode_fg(1:i34_fg(i),i)))/i34_fg(i)
    do j=1,ne
      xi=sum(xnd(elnode(1:i34(j),j)))/i34(j)
      yi=sum(ynd(elnode(1:i34(j),j)))/i34(j)
      dist=sqrt((xi_fg-xi)**2+(yi_fg-yi)**2) 
      if(dist<mdist) then
        mdist=dist
        l=j
      endif !dist
    enddo !j
    if(l==0) stop 'can not find mdist'
    ind2De(i)=l

    do k=kbe_fg(i),nvrt_fg
      ie=ind2De(i)
      zi_fg=0.0; isum=0
      do n=1,i34_fg(i)
        if(k/=nvrt_fg.and.zcor_fg(k,elnode_fg(n,i))==0.0) cycle
        zi_fg=zi_fg+zcor_fg(k,elnode_fg(n,i))
        isum=isum+1
      enddo
      zi_fg=zi_fg/isum

      l=0;  mdist=1.d6
      do m=kbe(ie),nvrt
        zi=0.0; isum=0
        do n=1,i34(ie)
          if(m/=nvrt.and.zcor(m,elnode(n,ie))==0.0) cycle
          zi=zi+zcor(m,elnode(n,ie))
          isum=isum+1
        enddo
        zi=zi/isum
        dist=abs(zi_fg-zi)
        if(dist<mdist) then
          mdist=dist
          l=m
        endif !dist
      enddo  !m
      if(l==0) stop 'can not find mdist in zcor'
      ind3De(k,i)=l
      if(k==kbe_fg(i)) ind3De(1:kbe_fg(i)-1,i)=l
    enddo !k
  enddo !i

  !read data from hotstart.in.0
  open(21,file='hotstart.in.0',form='unformatted',status='old')
 
  read(21)rtmp,itmp1,itmp2 
  do i=1,ne
    read(21)itmp1,itmp2,((tr_el(k,j,i),j=1,ntracers+1),k=1,nvrt)
  enddo
  do i=1,ns
    read(21)itmp,itmp,((rtmp,j=1,4),k=1,nvrt)
  enddo
  do i=1,np
    read(21)itmp,rtmp,itmp,((tr_nd(k,j,i),j=1,2*ntracers+7),k=1,nvrt)
  enddo
  !do i=1,ne
  !  read(21)itmp,(tr_sed(j,i),j=1,40)
  !enddo
  close(21)
 
  !write data for new hotstart 
  open(22,file='hotstart.in',form='unformatted',status='replace')
  write(22)0.d0,0,1
  do i=1,ne_fg
    write(22)i,0,((tr_el(ind3De(k,i),j,ind2De(i)),j=1,ntracers+1),k=1,nvrt_fg)
  enddo
  do i=1,ns_fg
    write(22)i,0,(0.d0,0.d0,0.d0,0.d0,j=1,nvrt_fg)
  enddo
  do i=1,np_fg
    write(22)i,0.d0,0,((tr_nd(ind3Dp(k,i),j,ind2Dp(i)),j=1,2*ntracers+7),k=1,nvrt_fg)
  enddo
  !do i=1,ne_fg
  !  write(22)i,(tr_sed(j,ind2De(i)),j=1,40)
  !enddo
  

  close(22)
  stop

end program Interp_hotstart


subroutine read_grid
!----------------------------------------------------------------
!read backgroud grid information 
!----------------------------------------------------------------
  use var_glbl,only : ne,np,ns,ivcor,nvrt,nxq,i34,elnode,nne,indel, &
                    & elside,ic3,isdel,isidenode,xnd,ynd,dp,kbp,kbe,zcor, &
                    & ztot,sigma,sigma_lcl,kz,h0,h_s,h_c,theta_b,theta_f,xctr,yctr

  !local variables
  implicit none
  integer, parameter :: mnei=20,mns=150000 !max. 
  integer :: i,j,k,l,ie,nd,nd1,nd2,idry,istat
  real(4) :: tmp

  !read bg grid info
  open(11,file='bg.gr3',status='old')
  read(11,*); read(11,*)ne,np

  open(12,file='vgrid.bg',status='old')
  read(12,*)ivcor; read(12,*)nvrt
  !if(ivcor/=1) stop 'ivcor/=1 in vgrid.bg'

  !alloc variables
  allocate(xnd(np),ynd(np),dp(np),i34(ne),elnode(4,ne),nne(np),indel(mnei,np), &
         & ic3(4,ne),elside(4,ne),isdel(2,mns),isidenode(2,mns),zcor(nvrt,np), &
         & kbp(np),kbe(ne),ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np), &
         & xctr(ne),yctr(ne),stat=istat)
  if(istat/=0) stop 'failed in alloc. xnd'
  
  do i=1,np
    read(11,*)j,xnd(i),ynd(i),dp(i)
  enddo
  do i=1,ne
    read(11,*)j,i34(i),(elnode(j,i),j=1,i34(i))
  enddo
  close(11); close(12);
  
  !read vgrid.bg
  call get_vgrid('vgrid.bg',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
  do i=1,np
    if(ivcor==2) then
      tmp=dp(i)
      call zcor_SZ(real(tmp,4),0.0,real(h0,4),h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,zcor(:,i),idry,kbp(i))
    elseif(ivcor==1) then
      zcor(kbp(i):nvrt,i)=dp(i)*sigma_lcl(kbp(i):nvrt,i)
    endif !ivcor
  enddo !i

  !for kbe
  xctr=0.0; yctr=0.0
  do i=1,ne
    kbe(i)=minval(kbp(elnode(1:i34(i),i)))
    do j=1,i34(i)
      xctr(i)=xctr(i)+xnd(elnode(j,i))/i34(i)
      yctr(i)=yctr(i)+ynd(elnode(j,i))/i34(i)
    enddo !j
  enddo !i

  !cyclic index
  do k=3,4 !elem. type
    do i=1,k  !local index
      do j=1,k-1 !offset
        nxq(j,i,k)=i+j
        if(nxq(j,i,k)>k) nxq(j,i,k)=nxq(j,i,k)-k
        if(nxq(j,i,k)<1.or.nxq(j,i,k)>k) then
          write(*,*)'INIT: nx wrong',i,j,k,nxq(j,i,k)
          stop
        endif
      enddo !j
    enddo !i
  enddo !k

  !for nne and indel
  nne=0
  do i=1,ne
    do j=1,i34(i)
      nd=elnode(j,i)
      nne(nd)=nne(nd)+1
      if(nne(nd)>mnei) stop 'nne(nd)>mnei in bg.gr3'
      indel(nne(nd),nd)=i
    enddo
  enddo !i

  !for ic3
  do i=1,ne
    do j=1,i34(i)
      ic3(j,i)=0 !index for bnd sides
      nd1=elnode(nxq(1,j,i34(i)),i)
      nd2=elnode(nxq(2,j,i34(i)),i)
      do k=1,nne(nd1)
        ie=indel(k,nd1)
        if(ie/=i.and.(elnode(1,ie)==nd2.or.elnode(2,ie)==nd2.or.elnode(3,ie)==nd2.or. &
          &  (i34(ie)==4.and.elnode(4,ie)==nd2))) ic3(j,i)=ie
      enddo 
    enddo!j
  enddo !i

  !ns,elside,isdel,isidenode
  ns=0
  do i=1,ne
    do j=1,i34(i)
      nd1=elnode(nxq(1,j,i34(i)),i)
      nd2=elnode(nxq(2,j,i34(i)),i)
      if(ic3(j,i)==0.or.i<ic3(j,i)) then !new side
        ns=ns+1
        if(ns>mns) stop 'ns>mns in bg.gr3'
        elside(j,i)=ns
        isdel(1,ns)=i
        isdel(2,ns)=ic3(j,i)
        isidenode(1,ns)=nd1
        isidenode(2,ns)=nd2
        if(ic3(j,i)/=0) then !old internal side
          ie=ic3(j,i)
          l=0
          do k=1,i34(ie)
            if(ic3(k,ie)==i) then
              l=k
              exit
            endif
          enddo !k
          if(l==0) stop 'wrong, can not find ic3(k,ie)==i' 
          elside(l,ie)=ns
        endif !ic3
      endif !ic3
    enddo !j
  enddo !i

  return
end subroutine read_grid

subroutine read_grid_fg
!----------------------------------------------------------------
!read fg grid information 
!----------------------------------------------------------------
  use var_glbl,only : ne_fg,np_fg,ns_fg,ivcor_fg,nvrt_fg,nxq,i34_fg,elnode_fg,nne_fg,indel_fg, &
                    & elside_fg,ic3_fg,isdel_fg,isidenode_fg,xnd_fg,ynd_fg,dp_fg,kbp_fg,kbe_fg,zcor_fg, &
                    & ztot_fg,sigma_fg,sigma_lcl_fg,kz_fg,h0,h_s_fg,h_c_fg,theta_b_fg,theta_f_fg,xctr_fg,yctr_fg

  !local variables
  implicit none
  integer, parameter :: mnei=20,mns=150000 !max. 
  integer :: i,j,k,l,ie,nd,nd1,nd2,idry,istat

  !read bg grid info
  open(11,file='fg.gr3',status='old')
  read(11,*); read(11,*)ne_fg,np_fg

  open(12,file='vgrid.fg',status='old')
  read(12,*)ivcor_fg; read(12,*)nvrt_fg

  !alloc variables
  allocate(xnd_fg(np_fg),ynd_fg(np_fg),dp_fg(np_fg),i34_fg(ne_fg),elnode_fg(4,ne_fg),nne_fg(np_fg),indel_fg(mnei,np_fg), &
         & ic3_fg(4,ne_fg),elside_fg(4,ne_fg),isdel_fg(2,mns),isidenode_fg(2,mns),zcor_fg(nvrt_fg,np_fg), &
         & kbp_fg(np_fg),kbe_fg(ne_fg),ztot_fg(nvrt_fg),sigma_fg(nvrt_fg),sigma_lcl_fg(nvrt_fg,np_fg), &
         & xctr_fg(ne_fg),yctr_fg(ne_fg),kz_fg,h_s_fg,h_c_fg,theta_b_fg,theta_f_fg,stat=istat)
  if(istat/=0) stop 'failed in alloc. xnd_fg'
  
  do i=1,np_fg
    read(11,*)j,xnd_fg(i),ynd_fg(i),dp_fg(i)
    dp_fg(i)=max(dp_fg(i),0.5)
  enddo
  do i=1,ne_fg
    read(11,*)j,i34_fg(i),(elnode_fg(j,i),j=1,i34_fg(i))
  enddo
  close(11); close(12);
  
  !read vgrid.bg
  call get_vgrid('vgrid.fg',np_fg,nvrt_fg,ivcor_fg,kz_fg,h_s_fg,h_c_fg,theta_b_fg,theta_f_fg,ztot_fg,sigma_fg,sigma_lcl_fg,kbp_fg)
  do i=1,np_fg
    if(ivcor_fg==2) then
      call zcor_SZ(real(dp_fg(i),4),0.0,real(h0,4),h_s_fg,h_c_fg,theta_b_fg,theta_f_fg,kz_fg,nvrt_fg,ztot_fg,sigma_fg,zcor_fg(:,i),idry,kbp_fg(i))
    elseif(ivcor_fg==1) then
      zcor_fg(kbp_fg(i):nvrt_fg,i)=dp_fg(i)*sigma_lcl_fg(kbp_fg(i):nvrt_fg,i)
    endif !ivcor
  enddo !i

  !for kbe
  xctr_fg=0.0; yctr_fg=0.0
  do i=1,ne_fg
    kbe_fg(i)=minval(kbp_fg(elnode_fg(1:i34_fg(i),i)))
    do j=1,i34_fg(i)
      xctr_fg(i)=xctr_fg(i)+xnd_fg(elnode_fg(j,i))/i34_fg(i)
      yctr_fg(i)=yctr_fg(i)+ynd_fg(elnode_fg(j,i))/i34_fg(i)
    enddo !j
  enddo !i

  !for nne and indel
  nne_fg=0
  do i=1,ne_fg
    do j=1,i34_fg(i)
      nd=elnode_fg(j,i)
      nne_fg(nd)=nne_fg(nd)+1
      if(nne_fg(nd)>mnei) stop 'nne(nd)>mnei in fg.gr3'
      indel_fg(nne_fg(nd),nd)=i
    enddo
  enddo !i

  !for ic3
  do i=1,ne_fg
    do j=1,i34_fg(i)
      ic3_fg(j,i)=0 !index for bnd sides
      nd1=elnode_fg(nxq(1,j,i34_fg(i)),i)
      nd2=elnode_fg(nxq(2,j,i34_fg(i)),i)
      do k=1,nne_fg(nd1)
        ie=indel_fg(k,nd1)
        if(ie/=i.and.(elnode_fg(1,ie)==nd2.or.elnode_fg(2,ie)==nd2.or.elnode_fg(3,ie)==nd2.or. &
          &  (i34_fg(ie)==4.and.elnode_fg(4,ie)==nd2))) ic3_fg(j,i)=ie
      enddo 
    enddo!j
  enddo !i

  !ns,elside,isdel,isidenode
  ns_fg=0
  do i=1,ne_fg
    do j=1,i34_fg(i)
      nd1=elnode_fg(nxq(1,j,i34_fg(i)),i)
      nd2=elnode_fg(nxq(2,j,i34_fg(i)),i)
      if(ic3_fg(j,i)==0.or.i<ic3_fg(j,i)) then !new side
        ns_fg=ns_fg+1
        if(ns_fg>mns) stop 'ns>mns in fg.gr3'
        elside_fg(j,i)=ns_fg
        isdel_fg(1,ns_fg)=i
        isdel_fg(2,ns_fg)=ic3_fg(j,i)
        isidenode_fg(1,ns_fg)=nd1
        isidenode_fg(2,ns_fg)=nd2
        if(ic3_fg(j,i)/=0) then !old internal side
          ie=ic3_fg(j,i)
          l=0
          do k=1,i34_fg(ie)
            if(ic3_fg(k,ie)==i) then
              l=k
              exit
            endif
          enddo !k
          if(l==0) stop 'wrong, can not find ic3(k,ie)==i' 
          elside_fg(l,ie)=ns_fg
        endif !ic3
      endif !ic3
    enddo !j
  enddo !i

  return
end subroutine read_grid_fg
