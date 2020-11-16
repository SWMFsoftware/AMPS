!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==========================================================================
module PT_wrapper

  use CON_coupler, ONLY: OH_,IH_,PT_,SC_, Couple_CC,Grid_C, iCompSourceCouple
  use CON_time
  use,intrinsic :: ieee_arithmetic

  implicit none

  private ! except

  ! Coupling with CON
  public:: PT_set_param
  public:: PT_init_session
  public:: PT_run
  public:: PT_save_restart
  public:: PT_finalize

  ! Point coupling
  public:: PT_get_grid_info 
  public:: PT_find_points

  ! GM coupling
  public:: PT_put_from_gm  

  ! OH coupling
  public:: PT_put_from_oh
  public:: PT_get_for_oh
  public:: PT_put_from_oh_dt

  ! IH coupling
  public:: PT_put_from_ih
  public:: PT_put_from_ih_dt

  ! SC coupling
  public:: PT_put_from_sc 
  public:: PT_put_from_sc_dt

  !codes describeing status of coupling with particular components of the SWMF 
  integer:: IhCouplingCode
  integer:: OhCouplingCode
  integer:: ScCouplingCode 


  !coupling operation counter (need for debugging)
  integer::nRecvFromOH=0
  integer::nSentToOH=0

contains

  subroutine PT_set_param(CompInfo, TypeAction)

    use CON_comp_info
    use ModReadParam

    integer :: iComm,iProc,nProc, nThread

    ! Arguments
    type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp
    character (len=*), intent(in)     :: TypeAction ! What to do

    ! Contains the PARAM.in segment
    character(len=lStringLine), allocatable :: StringLineF_I(:)

    character (len=*), parameter :: NameSub='PT_set_param'
    character (len=2) ComponentName

    !-------------------------------------------------------------------------
    ComponentName=CompInfo%name
    call amps_set_component_name(ComponentName)  

    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.true., &
            NameVersion='AMPS', &
            Version    =1.0)

    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc, &
            nThread=nThread)
       call AMPS_SetMpiCommunicator(iComm, iProc, nProc, nThread)

    case('CHECK')
       ! AMPS could check now the input parameters for consistency

    case('READ')
       ! get section of PARAM.in that contains the PT module
       allocate(StringLineF_I(i_line_read()+1:n_line_read()))
       call read_text(StringLineF_I)
       if(n_line_read()-i_line_read()+1 > 0) then
          call amps_read_param(StringLineF_I, n_line_read()-i_line_read(), &
               lStringLine,iProc)
       end if
    case('STDOUT')
       ! call AMPS_set_stdout

    case('FILEOUT')
       ! call AMPS_set_fileout

    case('GRID')
       ! Grid info depends on BATSRUS

    case default
       call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')
    end select

  end subroutine PT_set_param

  !============================================================================

  subroutine PT_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PT_init_session'
    integer::code

    Grid_C(PT_)%TypeCoord='HGI' 

    if (DoTimeAccurate) then 
      code=1
    else 
      code=0
    end if 
    
    call amps_init_session(iSession,TimeSimulation,code)
  end subroutine PT_init_session

  !============================================================================

  subroutine PT_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PT_finalize'
    !-------------------------------------------------------------------------
    call AMPS_finalize

  end subroutine PT_finalize

  !============================================================================

  subroutine PT_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PT_save_restart'
    !-------------------------------------------------------------------------
    !!! PT should save restart files !!!

  !  call amps_save_restart()
    
  end subroutine PT_save_restart

  !============================================================================

  subroutine PT_run(TimeSimulation, TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='PT_run'
    real:: xEarth(3)
    !-------------------------------------------------------------------------

    !update the location of the Earth in the coupled on the AMPS side when boupling with IH as active 
    if (IhCouplingCode==1) then 
      call GetEarthLocation(xEarth) 
      call set_earth_locaton_hgi(xEarth)
    end if

    !call AMPS  
    call AMPS_timestep(TimeSimulation, TimeSimulationLimit) 
  end subroutine PT_run

  !============================================================================

  subroutine PT_get_grid_info(nDimOut, iGridOut, iDecompOut)

    ! Provide information about AMPS grid. Set number of ion fluids.

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index (increases with AMR)
    integer, intent(out):: iDecompOut ! decomposition index

    integer:: nVarCouple
    integer:: nCommunicatedFluids
    
    character(len=*), parameter :: NameSub = 'PT_get_grid_info'
    !--------------------------------------------------------------------------
    if(iCompSourceCouple /= PT_)then
       nVarCouple = Grid_C(iCompSourceCouple)%nVar
       !write(*,*)'!!!', NameSub, i_proc(), ' nVarCouple=', nVarCouple

       IhCouplingCode=0
       OhCouplingCode=0
       ScCouplingCode=0

       nCommunicatedFluids=nVarCouple / 5 

       if (Couple_CC(OH_,PT_)%DoThis) then 
         OhCouplingCode=1 
         Grid_C(PT_)%TypeCoord='HGI'
       end if 


       if (Couple_CC(IH_,PT_)%DoThis) then 
         IhCouplingCode=1 
         Grid_C(PT_)%TypeCoord='HGI'
       end if


       if (Couple_CC(SC_,PT_)%DoThis) ScCouplingCode=1

       if (iCompSourceCouple == IH_) nCommunicatedFluids=1 
       if (iCompSourceCouple == SC_) nCommunicatedFluids=1 

       ! Pass number of fluids to this incorrectly named subroutine
       call amps_from_oh_init(nCommunicatedFluids,OhCouplingCode,IhCouplingCode)
    end if
    
    nDimOut    = 3
    iGridOut   = 1
    iDecompOut = 1

    call amps_mesh_id(iDecompOut)

  end subroutine PT_get_grid_info

  !============================================================================
  subroutine PT_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    integer, intent(in) :: nDimIn                ! dimension of position vector
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position

    ! Find array of points and return processor indexes owning them
    ! Could be generalized to return multiple processors...

    real:: Xyz_D(3) = 0.0
    integer:: iPoint, iProcFound

    character(len=*), parameter:: NameSub = 'PT_find_points'
    !--------------------------------------------------------------------------
    do iPoint = 1, nPoint
      Xyz_D(:)=0.0
      Xyz_D(1:nDimIn) = Xyz_DI(:,iPoint)
      
      call amps_get_point_thread_number(iProcFound, Xyz_D)
      iProc_I(iPoint) = iProcFound
      
    end do

  end subroutine PT_find_points

  !============================================================================

  subroutine PT_put_from_gm( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    real, intent(out), optional, allocatable:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter :: NameSub='PT_put_from_gm'
    !--------------------------------------------------------------------------
    if(present(Pos_DI))then
       ! set number of grid points on this processor
       call amps_get_center_point_number(nPoint)

       ! allocate position array
       allocate(Pos_DI(3,nPoint))

       ! get point positions from AMPS
       call amps_get_center_point_coordinates(Pos_DI) 

    elseif(present(Data_VI))then
       call amps_recieve_batsrus2amps_center_point_data(&
            NameVar//char(0), nVar, Data_VI, iPoint_I)
    else
       call CON_stop(NameSub//': neither Pos_DI nor Data_VI are present!')
    end if

  end subroutine PT_put_from_gm

  !============================================================================

  subroutine PT_put_from_oh( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data
    integer::i,j,jj  

    real, intent(out), optional, allocatable:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter :: NameSub='PT_put_from_oh'
    !--------------------------------------------------------------------------
    nRecvFromOH=nRecvFromOH+1 

    if(present(Pos_DI))then
       ! set number of grid points on this processor
       call amps_get_center_point_number(nPoint)
       
       ! allocate position array
       allocate(Pos_DI(3,nPoint))
       
       ! get point positions from AMPS
       call amps_get_center_point_coordinates(Pos_DI) 

    elseif(present(Data_VI))then

    do i = 1,nVar
      do j=1,nPoint
       jj=iPoint_I(j) 

       if (jj>0) then
         if (ieee_is_nan(Data_VI(i,jj))) then
           call CON_stop(NameSub//': nan')
         end if

         if (.not.ieee_is_finite(Data_VI(i,jj))) then
           call CON_stop(NameSub//': not finite')
         end if
        end if 
      end do
    end do



       call amps_recieve_batsrus2amps_center_point_data(&
            NameVar//char(0), nVar, Data_VI, iPoint_I)
    else
       call CON_stop(NameSub//': neither Pos_DI nor Data_VI are present!')
    end if

  end subroutine PT_put_from_oh

  subroutine PT_put_from_ih( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    real, intent(out), optional, allocatable:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter :: NameSub='PT_put_from_ih'
    !--------------------------------------------------------------------------
    if(present(Pos_DI))then
       ! set number of grid points on this processor
       call amps_get_center_point_number_ih(nPoint)

       ! allocate position array
       allocate(Pos_DI(3,nPoint))

       ! get point positions from AMPS
       call amps_get_center_point_coordinates_ih(Pos_DI)

    elseif(present(Data_VI))then
       call amps_recieve_batsrus2amps_center_point_data_ih(&
            NameVar//char(0), nVar, Data_VI, iPoint_I)
    else
       call CON_stop(NameSub//': neither Pos_DI nor Data_VI are present!')
    end if

  end subroutine PT_put_from_ih

  subroutine PT_put_from_sc(NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI) 
    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    real, intent(out), optional, allocatable:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter :: NameSub='PT_put_from_sc'
    !--------------------------------------------------------------------------
    if (present(Pos_DI)) then
       ! set number of grid points on this processor
       call amps_get_center_point_number_sc(nPoint)

       ! allocate position array
       allocate(Pos_DI(3,nPoint))

       ! get point positions from AMPS
       call amps_get_center_point_coordinates_sc(Pos_DI)
    elseif (present(Data_VI)) then
       call amps_recieve_batsrus2amps_center_point_data_sc(NameVar//char(0), nVar, Data_VI, iPoint_I) 
    else
       call CON_stop(NameSub//': neither Pos_DI nor Data_VI are present!')
    end if
  end subroutine PT_put_from_sc 


  !============================================================================
  subroutine PT_put_from_oh_dt(Dt)

    real,    intent(in):: Dt
    character(len=*), parameter :: NameSub='PT_put_from_oh_dt'
    !--------------------------------------------------------------------------
    call amps_impose_global_time_step(Dt)

  end subroutine PT_put_from_oh_dt

  subroutine PT_put_from_ih_dt(Dt)

    real,    intent(in):: Dt
    character(len=*), parameter :: NameSub='PT_put_from_ih_dt'
    !--------------------------------------------------------------------------
!   call amps_impose_global_time_step(Dt)

  end subroutine PT_put_from_ih_dt

  subroutine PT_put_from_sc_dt(Dt)

    real,    intent(in):: Dt
    character(len=*), parameter :: NameSub='PT_put_from_sc_dt'
    !--------------------------------------------------------------------------
!   call amps_impose_global_time_step(Dt)

  end subroutine PT_put_from_sc_dt

  !============================================================================
  subroutine PT_get_for_oh(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    ! Get data from PT to OH

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    integer::i,j

    character(len=*), parameter :: NameSub='PT_get_for_oh'
    !--------------------------------------------------------------------------
    nSentToOH=nSentToOH+1 

    call amps_send_batsrus2amps_center_point_data( &
         NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, Data_VI)

    do i = 1,nVarIn
      do j=1,nPoint  
       if (ieee_is_nan(Data_VI(i,j))) then   
         call CON_stop(NameSub//': nan')
       end if 

       if (.not.ieee_is_finite(Data_VI(i,j))) then
         call CON_stop(NameSub//': not finite')
       end if
      end do
    end do

  end subroutine PT_get_for_oh

end module PT_wrapper

subroutine ConvertX_HGI_HGR(xHGR,xHGI)

   use  CON_axes 

   real, intent(in) :: xHGI(3)  ! Position vectors
   real, intent(out) :: xHGR(3)  ! Position vectors

xHGR=matmul(HgrHgi_DD,xHGI)

!  xHGR=matmul(xHGI,HgiHgr_DD)
end subroutine ConvertX_HGI_HGR


subroutine ConvertX_HGR_HGI(xHGI,xHGR)
   use  CON_axes, ONLY:HgrHgi_DD

   implicit none

   real, intent(out) :: xHGI(3)  ! Position vectors
   real, intent(in) :: xHGR(3)  ! Position vectors



   xHGI=matmul(xHGR,HgrHgi_DD)
end subroutine ConvertX_HGR_HGI

subroutine ConvertVel_HGR_HGI(vHGI,vHGR,xHGR,TimeSim)
   use  CON_axes, ONLY:transform_velocity

   implicit none

   real, intent(out) :: vHGI(3) ! Position vector
  
   real, intent(in) :: xHGR(3)  ! Position vector
   real, intent(in) :: vHGR(3)  ! Velocity vector 
   real, intent(in) :: TimeSim  ! Simulation time

   vHGI=transform_velocity(TimeSim,vHGR,xHGR,'HGR','HGI')
end subroutine ConvertVel_HGR_HGI

subroutine GetEarthLocation(xEarthHgi) 
   use CON_axes
   use CON_time,ONLY: tSimulation

   implicit none

   real,intent(out)::xEarthHgi(3)

!   call set_hgi_gse_d_planet(tSimulation)
  
   ! Calculate the planet position in HGI                                                                                                   
   ! In GSE shifted to the center of the Sun the planet is at (-d,0,0)                                                                      
   xEarthHgi = matmul(HgiGse_DD, (/-cAU*SunEMBDistance, 0.0, 0.0/))
end subroutine GetEarthLocation





