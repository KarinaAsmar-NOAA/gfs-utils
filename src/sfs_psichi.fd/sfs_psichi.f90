!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM: SFS_PSICHI 
!   PRGMMR: ASMAR              ORG: EMC        DATE: 2024-12-18 
!
! ABSTRACT: THIS PROGRAM READS ISOBARIC WINDS ON GRID POINTS FROM 
! A GRIB2 FILE, TRANSFORMS THEM BACK TO SPECTURM SPACE USING THE
! NCEPLIBS-SP LIBRARY, AND COMPUTES PSI (STREAMFUNCTION), 
! CHI (VELOCITY POTENTIAL), VORTICITY AND DIVERGENCE. PSI AND CHI
! ARE SAVED AND OUTPUT IN GRIB2 FORMAT. IT IS BASED ON THE CFS_GENPSIANDCHI.FD
! UTILITY AUTHORED BY SAHA, CHUANG, AND MOORTHI.
!
! PROGRAM HISTORY LOG:
!   24-12-18   KARINA ASMAR           - ORIGINATOR
!
! USAGE:
!????????????????????????????????????????????????????????
!   INPUT FILES:  ??????????????????????
!     UNIT  11  GRIB FILE 
!     UNIT   5  READ *, IENST,IENSI (For ensemble message)
!
!   OUTPUT FILES:
!     UNIT  51  GRIB FILE  ???????????????
!
!   SUBPROGRAMS CALLED:
!     GETGB2 -- W3LIB ROUTINE ??????????????
!     PUTGB2 -- W3LIB ROUTINE ???????????????
!?????????????????????????????????????????????????????
! ATTRIBUTES:
!   LANGUAGE: FORTRAN
!
!$$$
      program sfs_psichi

      use grib_mod
      
      implicit none
!
      integer, parameter :: komax=100
      character*200 grbwnd,indwnd
      character*200 file,psichifile
      character*10 kyr,kmth
!
      real, allocatable ::  udata(:,:,:), vdata(:,:,:), vertlevs(:)
      real, allocatable ::  ui(:,:,:), vi(:,:,:), uo(:,:,:), vo(:,:,:)
      real, allocatable ::  div(:,:,:),  zo(:,:,:)
      real, allocatable ::  psio(:,:,:), so(:,:,:)
      real, allocatable ::  dummy1d(:), presslevs(:)
      character(len=1),allocatable,dimension(:) :: cgrib, cgrib2
!
      integer, parameter :: msk1=32000
!
      integer :: idim, jdim, kdim, ldim, jcap, idrt, lev, &
                 ierr, i, j, l, k, n, ij
      integer :: numfields,numlocal,maxlocal
      integer :: currlen=0,icount,itot,iseek
      integer :: lskip,lgrib,igdtlen,ipdtnum,ipdtlen
      integer :: idrtlen,idrtnum,ibmap,npt
      integer :: lengrib,listsec0(3),listsec1(13)
      integer :: listsec0_out(2),listsec1_out(13)
      integer, allocatable :: igdstmpl(:),ipdstmpl(:),idrtmpl(:)
      integer(4) :: max_bytes,igds(5),numcoord,coordlist
      
      type(gribfield) :: gfld

      logical*1,allocatable :: bmap(:)

       call getenv("PGBOUT",grbwnd)
       write(*,*) "grbwnd= ",grbwnd
!
       call getenv("PGIOUT",indwnd)
       write(*,*) "indwnd= ",indwnd
!
       call getenv("psichifile",psichifile)
       write(*,*) "psichifile= ",psichifile

       call baopenr(11,grbwnd,ierr)
         if(ierr.ne.0) then
         print *,'error opening file ',grbwnd
         call abort
       endif
       call baopenr(12,indwnd,ierr)
         if(ierr.ne.0) then
         print *,'error opening file ',indwnd
         call abort
       endif
!
       call baopenwt(51,psichifile,ierr)
         if(ierr.ne.0) then
         print *,'error opening file ',psichifile
         call abort
       endif

! allocate enough to get full dimensions
       allocate(udata(3000,3000,300))  
       allocate(vdata(3000,3000,300))
       allocate(vertlevs(300))

! Unpack GRIB2 as in WAFS/sorc/wafs_blending_0p25.fd/blending.f90
       icount=0
       iseek=0
       k=0
       do
          call skgb(11,iseek,msk1,lskip,lgrib)
          if (lgrib==0) exit    ! end loop at EOF or problem
          if (lgrib>currlen) then
             if (allocated(cgrib)) deallocate(cgrib)
             allocate(cgrib(lgrib),stat=ierr)
             currlen=lgrib
           endif
          call baread(11,lskip,lgrib,lengrib,cgrib)
             if (lgrib/=lengrib) then
             print *,' degrib2: IO Error.'
             call errexit(9)
          endif
          iseek=lskip+lgrib
          icount=icount+1

          call gb_info(cgrib,lengrib,listsec0,listsec1,  &
                       numfields,numlocal,maxlocal,ierr)
             if (ierr/=0) then
                 write(6,*) ' ERROR querying GRIB2 message = ',ierr
                 stop 10
             endif
             do n=1,numfields
             call gf_getfld(cgrib,lengrib,n,.true.,.true.,gfld,ierr)
                if (ierr/=0) then
                write(6,*) ' ERROR extracting field = ',ierr
                cycle
                endif
		
		! Get grid I,J dimensions and identifier for spectral truncation
  		idim=gfld%igdtmpl(8)
    		jdim=gfld%igdtmpl(9)
		!idrt=gfld%igdtnum  ! can also be gfld%griddef
		idrt=gfld%griddef

  		! get specs for output grib2 file
    		npt = gfld%ngrdpts
    		max_bytes = npt*4
		igdtlen = gfld%igdtlen
      		igdstmpl = gfld%igdtmpl
		ipdtnum = gfld%ipdtnum
  		ipdtlen	= gfld%ipdtlen
    		idrtnum = gfld%idrtnum
      		idrtlen = gfld%idrtlen
		idrtmpl = gfld%idrtmpl

 		! CATEGORY NO.: gfld%ipdtmpl(1) =2 METEO./MOMENTUM
   		! PARAMETER NO.: gfld%ipdtmpl(2) =2 UGRD =3 VGRD
     		! LEVEL ID: gfld%ipdtmpl(10) = 100 ISOBARIC SURFACE
       		! LEVEL VALUE: gfld%ipdtmpl(12) = MB

    		if ((gfld%ipdtmpl(1)==2) .and. (gfld%ipdtmpl(10)==100)) then
		  if (gfld%ipdtmpl(2)==2) then ! U/Wind
          	    k=k+1
	            vertlevs(k)=gfld%ipdtmpl(12)
	     	    print*,'vertical level',vertlevs(k)
	  	    do j=1,jdim
     		      do i=1,idim
	 		udata(i,j,k)=gfld%fld((j-1)*idim+i)
	  	      enddo
     		    enddo
	 	  endif
     		  if (gfld%ipdtmpl(2)==3) then ! V/Wind
	  	    do j=1,jdim
     		      do i=1,idim
	 		vdata(i,j,k)=gfld%fld((j-1)*idim+i)
	  	      enddo
     		    enddo
      		  endif
    	        endif
    	     enddo  ! end do numfields
!       enddo

 !       call gf_free(gfld)

 ! Allocate arrays, compute spectral max. wave for truncation
 ! and perform spectral truncation of winds

 ! Total isobaric levels 
       ldim=k
       allocate(presslevs(ldim))
       do l=1,ldim
         presslevs(l)=vertlevs(l)
	 print*,'final press mb',l,presslevs(l)
       enddo

       allocate (ui(idim,jdim,ldim))
       allocate (vi(idim,jdim,ldim))
       allocate (uo(idim,jdim,ldim))
       allocate (vo(idim,jdim,ldim))
       allocate (div(idim,jdim,ldim))
       allocate (zo(idim,jdim,ldim))
       allocate (psio(idim,jdim,ldim))
       allocate (so(idim,jdim,ldim))
 	do l=1,ldim
  	  do j=1,jdim
   	    do i=1,idim
    	      ui(i,j,l)=udata(i,j,l)
              vi(i,j,l)=vdata(i,j,l)
    	    enddo
   	  enddo
  	enddo

	print*,'udata',udata(idim/2,jdim/2,ldim/2)
	print*,'vdata',vdata(idim/2,jdim/2,ldim/2)

	deallocate(udata)
 	deallocate(vdata)

	!idrt=0    !!! check how to get this right.... !=0 yields psi=0
    	if(idrt == 0)then
	   jcap = (jdim-3)/2
	 else
	   jcap = jdim-1
	 end if

	print*,'ui',ui(idim/2,jdim/2,ldim/2)
	print*,'vi',vi(idim/2,jdim/2,ldim/2)
	print*,'jcap',jcap
	print*,'idrt',idrt
        call sptrunv(0,jcap,idrt,idim,jdim,idrt,idim,jdim,ldim,          &
                     0,0,0,0,0,0,0,0,ui(1,1,1),vi(1,1,1),		 &
                    .false.,uo(1,1,1),vo(1,1,1),.false.,div,zo,.true.	 &
                    ,psio(1,1,1),so(1,1,1))

	deallocate(ui)
 	deallocate(vi)
  	deallocate(uo)
   	deallocate(vo)
    	deallocate(div)
     	deallocate(zo)
        do l=1,ldim
	  print*,'psio',l,psio(idim/2,jdim/2,l)
   	  print*,'so',l,so(idim/2,jdim/2,l)
	enddo

! Write GRIB2 file for each mb level
      do l=1,ldim
        lev=presslevs(l)
        allocate(cgrib2(max_bytes))

!==>initialize new GRIB2 message and pack
! GRIB2 sections 0 (Indicator Section) and 1 (Identification Section)
        listsec0_out(1)=0 ! Discipline-GRIB Master Table Number (see Code Table 0.0)
        listsec0_out(2)=2 ! GRIB Edition Number (currently 2)
    	listsec1_out(1)=7 ! Id of orginating centre (Common Code Table C-1)  !!! CHECK**********
    	listsec1_out(2)=4 !"EMC"! Id of orginating sub-centre (local table)/Table C of ON388 CHECK *************
    	listsec1_out(4)=1    ! per Brent! GRIB Local Tables Version Number (Code Table 1.1)  **********
    	listsec1_out(5)=1    ! Significance of Reference Time (Code Table 1.2)  *********************
    	listsec1_out(10) = 0 ! Reference Time - Minute
    	listsec1_out(11) = 0 ! Reference Time - Second
    	listsec1_out(12) = 0 ! Production status of data (Code Table 1.3)
    	listsec1_out(13) = 1 ! Type of processed data (Code Table 1.4)
                     ! 0 for analysis products and 1 for forecast products
		     
        call gribcreate(cgrib2,max_bytes,listsec0_out,listsec1_out,ierr)

        igds(1)=0   !Source of grid definition (see Code Table 3.0)
    	igds(2)=npt !Number of grid points in the defined grid.
    	igds(3)=0   !Number of octets needed for each additional grid points definition
    	igds(4)=0   !Interpretation of list for optional points definition (Code Table 3.11)
    	igds(5)=idrt !Grid Definition Template Number (Code Table 3.1)
    
	call addgrid(cgrib2,max_bytes,igds,igdstmpl,igdtlen,0,1,ierr)

! ADD STREAMFUNCTION
     	allocate(ipdstmpl((ipdtlen)))
        ipdstmpl(1) = 2    ! ==> parameter category (see Code Table 4.1)
	ipdstmpl(2) = 4    ! ==> parameter number (see Code Tavle 4.2)
	ipdstmpl(3) = 2    ! ==> type of generating process (see Code Table 4.3)
	ipdstmpl(4) = 0    ! ==> background generating process identifier 
	ipdstmpl(5) = 96   ! ==> analysis or forecast generating process identifier (defined by originating Center)
 	ipdstmpl(6) = 0    ! ==> hours of observational data cutoff after reference time
  	ipdstmpl(7) = 0    ! ==> minutes of observational data cutoff after reference time
	ipdstmpl(8) = 1    ! ==> indicator of unit of time range (see Code Table 4.4)
			   ! ==>ifield4(9):forecast time in units defined by ifield4(8) 
	ipdstmpl(10) = 100 ! ==> type of first fixed surface (see Code Table 4.5)
    	ipdstmpl(11) = 0   ! ==> scale factor of first fixed surface
	ipdstmpl(12) = lev ! ==> scaled value of first fixed surface
 	ipdstmpl(13) = 255 ! ==> type of second fixed surface(See Code Table 4.5)
	ipdstmpl(14) = 0   ! ==> scale factor of second fixed surface
	ipdstmpl(15) = 0   ! == > scaled value of second fixed surface
    	if(size(ipdstmpl)>=16) then
!==> ifield4(16):Statistical process used within the spatial area (see Code Table 4.10)
	  ipdstmpl(17) = 3 ! ==> Type of spatial processing
          ipdstmpl(18) = 1 ! ==> Number of data points used in spatial processing
    	end if

	numcoord = 0
        coordlist = 0.
	ibmap = 255
	
	allocate(bmap(npt))

    	print *, "npt=",npt
    	print *, "ipdtnum=",ipdtnum,ipdtlen,ipdstmpl
    	print *, "coordlist=",coordlist,numcoord
    	print *, "idrtnum=",idrtnum,idrtlen,idrtmpl
	print*,  "ipdstmpl",ipdstmpl
	print*,  "psio", psio(idim/2,jdim/2,ldim/2)

 	! build the array field for grib2
  	allocate(dummy1d(idim*jdim))
   	ij=1
        do j=1,jdim
	  do i=1,idim
     	    dummy1d(ij)=psio(i,j,l)
	    ij=ij+1
     	  enddo
	enddo

      	print*,'dummy1d',shape(dummy1d)
       	print*,'dimensions',idim,jdim,ldim
	print*,'grid points',idim*jdim,npt
        
        call addfield(cgrib2,max_bytes,ipdtnum,ipdstmpl,ipdtlen,  &
                      coordlist,numcoord,idrtnum,idrtmpl,	  &
                      idrtlen,dummy1d,npt,ibmap,bmap,ierr)
 	print*,'add field err',ierr
  	
   	call gribend(cgrib2,max_bytes,lengrib,ierr)
	print*,'gribend err',ierr
    
        call wryte(51, lengrib, cgrib2)

 	deallocate(cgrib2)
       enddo   ! end do for ldim grib2
      enddo    ! enddo for everything???
      stop
      end
