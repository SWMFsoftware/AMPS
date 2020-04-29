!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==========================================================================
module PT_wrapper

  use CON_coupler, ONLY: PT_, Grid_C, iCompSourceCouple

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
    !-------------------------------------------------------------------------
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
    
  end subroutine PT_save_restart

  !============================================================================

  subroutine PT_run(TimeSimulation, TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='PT_run'
    !-------------------------------------------------------------------------
    call AMPS_timestep(TimeSimulation, TimeSimulationLimit) 

  end subroutine PT_run

  !============================================================================

  subroutine PT_get_grid_info(nDimOut, iGridOut, iDecompOut)

    ! Provide information about AMPS grid. Set number of ion fluids.

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index (increases with AMR)
    integer, intent(out):: iDecompOut ! decomposition index

    integer:: nVarCouple
    
    character(len=*), parameter :: NameSub = 'PT_get_grid_info'
    !--------------------------------------------------------------------------
    if(iCompSourceCouple /= PT_)then
       nVarCouple = Grid_C(iCompSourceCouple)%nVar
       !write(*,*)'!!!', NameSub, i_proc(), ' nVarCouple=', nVarCouple
       ! Pass number of fluids to this incorrectly named subroutine
       call amps_from_oh_init(nVarCouple / 5)
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

    real, intent(out), optional, allocatable:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter :: NameSub='PT_put_from_oh'
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

  end subroutine PT_put_from_oh
  !============================================================================
  subroutine PT_put_from_oh_dt(Dt)

    real,    intent(in):: Dt
    character(len=*), parameter :: NameSub='PT_put_from_oh_dt'
    !--------------------------------------------------------------------------
    call amps_impose_global_time_step(Dt)

  end subroutine PT_put_from_oh_dt

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

    character(len=*), parameter :: NameSub='PT_get_for_oh'
    !--------------------------------------------------------------------------
    ! Overly long ame that violates Fortran 90 rules !!!
    call amps_send_batsrus2amps_center_point_data( &
         NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, Data_VI)

  end subroutine PT_get_for_oh

end module PT_wrapper
