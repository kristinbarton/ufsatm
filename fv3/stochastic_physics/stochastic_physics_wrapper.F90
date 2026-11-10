!> @file
!> @brief Wrapper for stochastic physics in the FV3 model.
!> @author Dominikus Heinzeller @date 8/6/21
module stochastic_physics_wrapper_mod

  use machine, only: kind_phys

  implicit none

  ! For stochastic physics pattern generation
  !> grid latitude in radians, default to pi/2 to
  !> -pi/2 range, otherwise adj in subr is called
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlat
  !> grid longitude in radians, ok for both 0->2pi
  !> or -pi -> +pi ranges
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlon
  !> Array of SPPT (Stochastically Perturbed Parameterization Tendencies) weights
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sppt_wts
  !> Allocatable 3D array to store stochastic humidity weights
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: shum_wts
  !> Array of stochastic u wind components of the kinetic energy backscatter scheme.
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebu_wts
  !> Array of stochastic v wind components of the kinetic energy backscatter scheme.
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebv_wts
  !> Surface weights array used in stochastic physics calculations.
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sfc_wts
  !> Array to store stochastic physics perturbation weights.
  real(kind=kind_phys), dimension(:,:,:,:), allocatable, save :: spp_wts

  logical, save :: is_initialized = .false. !< logical if initialized

  integer, save :: lsoil = -999 !< number of soil layers

  !> total soil moisture
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: smc
  !> soil temperature
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: stc
  !> liquid soil moisture
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: slc
  !> vegetation fraction
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: vfrac
  !> maximum snow albedo, in fraciton
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: snoalb
  !> mean nir albedo with strong cosz dependency
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnsf
  !> mean nir albedo with weak cosz dependency
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnwf
  !> surface lw emissivity, in fraction
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: semis
  !> surface roughness length, in cm
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: zorll

  !real(kind=kind_phys), dimension(:,:),   allocatable, save :: stype
  integer, dimension(:,:),   allocatable, save :: stype !< soil type

  ! For cellular automata
  !> Sea surface temperature (SST) array, allocated dynamically.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: sst
  !> Land mask array for stochastic physics processes.
  !> The array is allocatable and saved across subroutine calls.  
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lmsk
  !> A 2D real array to store lake data.
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lake
  !> Allocatable array to store the u-component of wind (zonal wind) 
  !> with dimensions (latitude, longitude, vertical levels).
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: uwind
  !> Allocatable array to store the v-component of wind (velocity in the y-direction).
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: vwind
  !> Array to store the height values with physical kind precision.
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: height
  !> 2D real array to store the stochastic perturbations.
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: dx
  !> A 2D real array used to store 2-dim grid spacing variable
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: condition
  !> Allocatable 2D array to store deep convection coupling coefficients
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca_deep_cpl
  !> Array for storing turbulence coupling coefficients
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca_turb_cpl
  !> Array for shallow convection coupling coefficients
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca_shal_cpl
  !> Allocatable 2D array for coupling data, saved across subroutine calls
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca1_cpl
  !> 2D array for coupling data in stochastic physics module
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca2_cpl
  !> Allocatable 2D array for coupling data in stochastic physics
  !> The array is allocatable and saved across subroutine calls.
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca3_cpl


!----------------
! Public Entities
!----------------
! functions
  public stochastic_physics_wrapper

  contains

!-------------------------------
!  CCPP step
!-------------------------------
   !> Wrapper subroutine for stochastic physics in the FV3 model.
   !>
   !> This subroutine serves as a wrapper for the stochastic physics 
   !> processes in the FV3 atmospheric model. It interfaces with various 
   !> components of the model to apply stochastic perturbations.
   !>
   !> @param[in]  GFS_Control   Control structure for the GFS model.
   !> @param[in]  GFS_Statein   Input state variables for the GFS model.
   !> @param[in]  GFS_Grid      Grid structure for the GFS model.
   !> @param[inout]  GFS_Sfcprop   Surface properties for the GFS model.
   !> @param[inout]  GFS_Radtend   Radiative tendencies for the GFS model.
   !> @param[inout]  GFS_Coupling  Coupling structure for the GFS model.
   !> @param[inout]  Atm_block     Atmospheric block data.
   !> @param[out] ierr          Error return code.
   !>
   !> @author Dominikus Heinzeller @date 8/6/21
  subroutine stochastic_physics_wrapper (GFS_Control, GFS_Statein, GFS_Grid, GFS_Sfcprop, GFS_Radtend, GFS_Coupling, Atm_block, ierr)

#ifdef _OPENMP
    use omp_lib
#endif

    use GFS_typedefs,       only: GFS_control_type, GFS_statein_type, GFS_grid_type, GFS_sfcprop_type, GFS_radtend_type, GFS_coupling_type
    use mpp_mod,            only: FATAL, mpp_error
    use block_control_mod,  only: block_control_type
    use atmosphere_mod,     only: Atm, mygrid

    use stochastic_physics,           only: init_stochastic_physics, run_stochastic_physics
!    use cellular_automata_global_mod, only: cellular_automata_global
    use cellular_automata_sgs_mod,    only: cellular_automata_sgs
    use lndp_apply_perts_mod,         only: lndp_apply_perts
    !! use new cellular automata 
    use cellular_automata_mod
    use spectral_transforms, only: stochy_la2ga
    use fv_mp_mod, only: mp_reduce_min
    use kinddef, only: kind_io8

    implicit none

    type(GFS_control_type),   intent(inout) :: GFS_Control
    type(GFS_statein_type),   intent(in)    :: GFS_Statein
    type(GFS_grid_type),      intent(in)    :: GFS_Grid
    type(GFS_sfcprop_type),   intent(inout) :: GFS_Sfcprop
    type(GFS_radtend_type),   intent(inout) :: GFS_Radtend
    type(GFS_coupling_type),  intent(inout) :: GFS_Coupling
    type(block_control_type), intent(inout) :: Atm_block
    integer,                  intent(out)   :: ierr

    integer :: nthreads, nb, levs, maxblk, nblks, n, v, ixs, ixe
    logical :: param_update_flag

    !! cellular automata object, save makes sure it persists
    type(ca_grid_t), save :: CA_grid
    real(kind=kind_phys) :: grid_size
    real(kind=kind_io8), allocatable, save :: fv3_ca(:,:,:)

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif
    ierr = 0

    levs   = GFS_Control%levs
    maxblk = maxval(GFS_Control%blksz)
    nblks  = Atm_block%nblks

    ! Initialize

    initalize_stochastic_physics: if (.not. is_initialized) then

      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type > 0) .OR. GFS_Control%do_spp) then
         allocate(xlat(1:nblks,maxblk))
         allocate(xlon(1:nblks,maxblk))
         call transfer_field_to_stochastics(GFS_Control%blksz, GFS_Grid%xlat, xlat)
         call transfer_field_to_stochastics(GFS_Control%blksz, GFS_Grid%xlon, xlon)
        ! Initialize stochastic physics
        call init_stochastic_physics(levs, GFS_Control%blksz, GFS_Control%dtp, GFS_Control%sppt_amp,                                  &
            GFS_Control%input_nml_file, GFS_Control%fn_nml, GFS_Control%nlunit, xlon, xlat, GFS_Control%do_sppt, GFS_Control%do_shum, &
            GFS_Control%do_skeb, GFS_Control%lndp_type, GFS_Control%n_var_lndp, GFS_Control%use_zmtnblck, GFS_Control%skeb_npass,     &
            GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,    &
            GFS_Control%n_var_spp, GFS_Control%spp_var_list, GFS_Control%spp_prt_list, GFS_Control%spp_stddev_cutoff, GFS_Control%do_spp,                            &
            GFS_Control%ak, GFS_Control%bk, nthreads, GFS_Control%master, GFS_Control%communicator, ierr)
            if (ierr/=0)  then
                    write(6,*) 'call to init_stochastic_physics failed'
                    return
            endif
      end if
      if (GFS_Control%do_sppt) then
         allocate(sppt_wts(1:nblks,maxblk,1:levs))
      end if
      if (GFS_Control%do_shum) then
         allocate(shum_wts(1:nblks,maxblk,1:levs))
      end if
      if (GFS_Control%do_skeb) then
         allocate(skebu_wts(1:nblks,maxblk,1:levs))
         allocate(skebv_wts(1:nblks,maxblk,1:levs))
      end if
      if ( GFS_Control%do_spp ) then
         allocate(spp_wts(1:nblks,maxblk,1:levs,1:GFS_Control%n_var_spp))
         do n=1,GFS_Control%n_var_spp
           select case (trim(GFS_Control%spp_var_list(n)))
           case('pbl')
             GFS_Control%spp_pbl = 1
           case('sfc')
             GFS_Control%spp_sfc = 1
           case('mp')
             GFS_Control%spp_mp = 7
           case('rad')
             GFS_Control%spp_rad = 1
           case('gwd')
             GFS_Control%spp_gwd = 1
           case('cu_deep')
             GFS_Control%spp_cu_deep = 1
           end select
         end do
      end if
      if ( GFS_Control%lndp_type == 2 ) then
          allocate(sfc_wts(1:nblks,maxblk,1:GFS_Control%n_var_lndp))
          if ( (GFS_Control%lsm == GFS_Control%lsm_noah) .or. (GFS_Control%lsm == GFS_Control%lsm_noahmp)) then
            lsoil = GFS_Control%lsoil
          elseif (GFS_Control%lsm == GFS_Control%lsm_ruc) then
            lsoil = GFS_Control%lsoil_lsm
          endif
          allocate(smc   (1:nblks, maxblk, lsoil))
          do v = 1,GFS_Control%n_var_lndp
            select case (trim(GFS_Control%lndp_var_list(v)))
            case('smc')
              allocate(slc   (1:nblks, maxblk, lsoil))
              allocate(stype (1:nblks, maxblk))
            case('stc')
              allocate(stc   (1:nblks, maxblk, lsoil))
            case('vgf') 
              allocate(vfrac (1:nblks, maxblk))
            case('alb') 
              allocate(alnsf (1:nblks, maxblk))
              allocate(alnwf (1:nblks, maxblk))
            case('sal') 
              allocate(snoalb(1:nblks, maxblk))
            case('emi') 
              allocate(semis (1:nblks, maxblk))
            case('zol') 
              allocate(zorll (1:nblks, maxblk))
            endselect 
          enddo 
      endif


      if ( GFS_Control%lndp_type == 1 ) then ! this scheme sets perts once
         allocate(sfc_wts(1:nblks, maxblk, GFS_Control%n_var_lndp))
         call run_stochastic_physics(levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz,       &
                                     sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts,         &
                                     skebv_wts=skebv_wts, sfc_wts=sfc_wts,                              &
                                     spp_wts=spp_wts, nthreads=nthreads)
         ! Copy contiguous data back
         do nb=1,nblks
            GFS_Coupling%sfc_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
         end do
         deallocate(sfc_wts)
      end if
      ! Consistency check for cellular automata
      if(GFS_Control%do_ca)then
        if(GFS_Control%ca_sgs)then
           allocate(sst         (1:nblks, maxblk))
           allocate(lmsk        (1:nblks, maxblk))
           allocate(lake        (1:nblks, maxblk))
           allocate(uwind       (1:nblks, maxblk, 1:levs))
           allocate(vwind       (1:nblks, maxblk, 1:levs))
           allocate(height      (1:nblks, maxblk, 1:levs))
           allocate(condition   (1:nblks, maxblk))
           allocate(dx          (1:nblks, maxblk))
           allocate(ca_deep_cpl (1:nblks, maxblk))
           allocate(ca_turb_cpl (1:nblks, maxblk))
           allocate(ca_shal_cpl (1:nblks, maxblk))
        endif
        if(GFS_Control%ca_global)then
          ! Allocate contiguous arrays; no need to copy in (intent out)
          allocate(ca1_cpl (1:nblks, maxblk))
          allocate(ca2_cpl (1:nblks, maxblk))
          allocate(ca3_cpl (1:nblks, maxblk))
          allocate(fv3_ca(1:nblks, maxblk, GFS_Control%nca_g))
          fv3_ca = 0.0
        endif
      endif

      !! cellular automata init or restart -- where in init should this go?
      if (GFS_Control%do_ca .and. GFS_Control%ca_global .and. GFS_Control%nca_g > 0) then
        grid_size = minval(GFS_Grid%dx)
        call mp_reduce_min(grid_size)
        call initialize_cellular_automata(grid_size,            &
                                          GFS_Control%ncells_g, &
                                          GFS_Control%nca_g,    &
                                          CA_grid)
      endif

      is_initialized = .true.

    else initalize_stochastic_physics
      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type == 2) .OR. GFS_Control%do_spp) then
         call run_stochastic_physics(levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz, &
                                 sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
                                 spp_wts=spp_wts, nthreads=nthreads)
         ! Copy contiguous data back
         if (GFS_Control%do_sppt) then
            do nb=1,nblks
                GFS_Coupling%sppt_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = sppt_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_shum) then
            do nb=1,nblks
                GFS_Coupling%shum_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = shum_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_skeb) then
            do nb=1,nblks
                GFS_Coupling%skebu_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = skebu_wts(nb,1:GFS_Control%blksz(nb),:)
                GFS_Coupling%skebv_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = skebv_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_spp) then
            do n=1,GFS_Control%n_var_spp
               select case (trim(GFS_Control%spp_var_list(n)))
               case('pbl')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_pbl(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('sfc')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_sfc(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('mp')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_mp(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('gwd')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_gwd(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('rad')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_rad(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('cu_deep')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_cu_deep(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               end select
            end do
         end if

         if (GFS_Control%lndp_type == 2) then ! save wts, and apply lndp scheme
             do nb=1,nblks
                GFS_Coupling%sfc_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
             end do
 
             do nb=1,nblks
                ixs = GFS_control%chunk_begin(nb)
                ixe = GFS_control%chunk_end(nb)
                do v = 1,GFS_Control%n_var_lndp
                  ! used to identify locations with land model (=soil) 
                  if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                     smc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%smois(ixs:ixe,1:lsoil)
                  else  ! noah or noah-MP
                     smc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%smc(ixs:ixe,1:lsoil)
                  endif

                  select case (trim(GFS_Control%lndp_var_list(v)))
                  ! DH* is this correct? shouldn't this be slc ?
                  case('smc')
                      ! stype used to fetch soil params
                      stype(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%stype(ixs:ixe)
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                         slc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%sh2o(ixs:ixe,1:lsoil)
                      else  ! noah or noah-MP
                         slc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%slc(ixs:ixe,1:lsoil)
                      endif
                  case('stc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                         stc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%tslb(ixs:ixe,1:lsoil)
                      else ! noah or noah-MP 
                         stc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%stc(ixs:ixe,1:lsoil)
                      endif 
                  case('vgf')
                      if ( (GFS_Control%lsm == GFS_Control%lsm_noahmp) ) then 
                         ! assumes iopt_dveg = 4 (will be checked later)
                         vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%shdmax(ixs:ixe)
                      else ! ruc or noah-MP
                         vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%vfrac(ixs:ixe)
                      endif
                  case('alb')
                      alnsf(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%alnsf(ixs:ixe)
                      alnwf(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%alnwf(ixs:ixe)
                  case('sal')
                      snoalb(nb,1:GFS_Control%blksz(nb)) = GFS_Sfcprop%snoalb(ixs:ixe)
                  case('emi')
                      semis(nb,1:GFS_Control%blksz(nb))  = GFS_Radtend%semis(ixs:ixe)
                  case('zol')
                      zorll(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%zorll(ixs:ixe)
                  endselect
              enddo
             enddo


             param_update_flag = .false.
             ! noah and noah-MP treated differently, as global cycle doesn't overwrite shdmax for Noah-MP 
             ! determine whether land paramaters have been over-written to
             ! trigger applying perturbations (logic copied from GFS_driver),
             if ( (GFS_Control%lsm == GFS_Control%lsm_noah) .and. GFS_Control%nscyc >  0) then
                 if (mod(GFS_Control%kdt,GFS_Control%nscyc) == 1 ) then
                   param_update_flag = .true.
                 endif
             endif
             if ( ( GFS_Control%nscyc ==  0 .or. GFS_Control%lsm == GFS_Control%lsm_noahmp) .and. GFS_Control%first_time_step ) then
             ! call once at start of the forecast.
                    param_update_flag = .true.
             endif
              
             call lndp_apply_perts(GFS_Control%blksz, GFS_Control%lsm, GFS_Control%lsm_noah, GFS_Control%lsm_ruc,             &
                               GFS_Control%lsm_noahmp, GFS_Control%iopt_dveg, lsoil, GFS_Control%dtp, GFS_Control%kdt,        &
                               GFS_Control%n_var_lndp, GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,                  &
                               sfc_wts, xlon, xlat, stype, GFS_Control%pores, GFS_Control%resid,param_update_flag,            &
                               smc, slc, stc, vfrac, alnsf, alnwf, snoalb, semis, zorll, ierr)

             if (ierr/=0)  then
                    write(6,*) 'call to GFS_apply_lndp failed'
                    return
             endif

             do nb=1,nblks
                ixs = GFS_control%chunk_begin(nb)
                ixe = GFS_control%chunk_end(nb)

                do v = 1,GFS_Control%n_var_lndp

                  select case (trim(GFS_Control%lndp_var_list(v)))
                  case('smc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                           GFS_Sfcprop%smois(ixs:ixe,1:lsoil) = smc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                           GFS_Sfcprop%sh2o(ixs:ixe,1:lsoil)  = slc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      else  ! noah or noah-MP
                           GFS_Sfcprop%smc(ixs:ixe,1:lsoil) = smc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                           GFS_Sfcprop%slc(ixs:ixe,1:lsoil) = slc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      endif
                  case('stc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                           GFS_Sfcprop%tslb(ixs:ixe,1:lsoil)  = stc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      else ! noah or noah-MP 
                           GFS_Sfcprop%stc(ixs:ixe,1:lsoil) = stc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      endif 
                  case('vgf')
                      if ( (GFS_Control%lsm == GFS_Control%lsm_noahmp) ) then 
                        GFS_Sfcprop%shdmax(ixs:ixe)  = vfrac(nb,1:GFS_Control%blksz(nb))
                      else 
                        GFS_Sfcprop%vfrac(ixs:ixe)  = vfrac(nb,1:GFS_Control%blksz(nb))
                      endif
                  case('alb')
                      GFS_Sfcprop%alnsf(ixs:ixe)  = alnsf(nb,1:GFS_Control%blksz(nb))
                      GFS_Sfcprop%alnwf(ixs:ixe)  = alnwf(nb,1:GFS_Control%blksz(nb))
                  case('sal')
                      GFS_Sfcprop%snoalb(ixs:ixe) = snoalb(nb,1:GFS_Control%blksz(nb))
                  case('emi')
                      GFS_Radtend%semis(ixs:ixe)  = semis(nb,1:GFS_Control%blksz(nb))
                  case('zol')
                      GFS_Sfcprop%zorll(ixs:ixe)  = zorll(nb,1:GFS_Control%blksz(nb))
                  end select
                enddo 
            enddo
         endif ! lndp block
      endif ! if do* block

      if (GFS_Control%do_ca) then

       if(GFS_Control%ca_sgs)then
         do nb=1,nblks
             ixs = GFS_control%chunk_begin(nb)
             ixe = GFS_control%chunk_end(nb)
             sst        (nb,1:GFS_Control%blksz(nb)) = GFS_Sfcprop%tsfco(ixs:ixe)
             lmsk       (nb,1:GFS_Control%blksz(nb)) = GFS_Sfcprop%slmsk(ixs:ixe)
             lake       (nb,1:GFS_Control%blksz(nb)) = GFS_Sfcprop%lakefrac(ixs:ixe)
             condition  (nb,1:GFS_Control%blksz(nb)) = GFS_Coupling%condition(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb))
             ca_deep_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Coupling%ca_deep(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb))
             ca_turb_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Coupling%ca_turb(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb))
             ca_shal_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Coupling%ca_shal(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb))
          enddo
          call transfer_field_to_stochastics_3d(GFS_Control%blksz, GFS_Statein%ugrs, uwind)
          call transfer_field_to_stochastics_3d(GFS_Control%blksz, GFS_Statein%vgrs, vwind)
          call transfer_field_to_stochastics_3d(GFS_Control%blksz, GFS_Statein%phil, height)
         call transfer_field_to_stochastics(GFS_Control%blksz, GFS_Grid%dx, dx)
         !! call the update -- currently only global flag available
         !call update_cellular_automata(
         !       GFS_Control%kdt,            &
         !       GFS_Control%iseed_ca,       &
         !       GFS_Control%nfracseed,      &
         !       GFS_Control%nseed,          &
         !       GFS_Control%nspinup,        &
         !       GFS_Control%nlives,         &
         !       CA_grid
         !     )
         call cellular_automata_sgs(GFS_Control%kdt,GFS_control%dtp,GFS_control%restart,GFS_Control%first_time_step,              &
            sst,lmsk,lake,uwind,vwind,height,dx,condition,ca_deep_cpl,ca_turb_cpl,ca_shal_cpl, Atm(mygrid)%domain_for_coupler,nblks,      &
            Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy, levs,                        &
            GFS_Control%nthresh,GFS_Control%tile_num,GFS_Control%nca,GFS_Control%ncells,GFS_Control%nlives,                       &
            GFS_Control%nfracseed, GFS_Control%nseed,GFS_Control%iseed_ca,GFS_Control%ca_advect,                                  &
            GFS_Control%nspinup,GFS_Control%ca_trigger,Atm_block%blksz(1),GFS_Control%master,GFS_Control%communicator)
         ! Copy contiguous data back as needed
         do nb=1,nblks
             GFS_Coupling%ca_deep(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca_deep_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Coupling%ca_turb(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca_turb_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Coupling%ca_shal(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca_shal_cpl (nb,1:GFS_Control%blksz(nb))
         enddo
       endif
       if(GFS_Control%ca_global)then
          if (GFS_Control%nca_g > 0) then
              call update_cellular_automata(GFS_Control%kdt,       &
                                            GFS_Control%iseed_ca,  &
                                            GFS_Control%nfracseed, &
                                            GFS_Control%nseed_g,   &
                                            GFS_Control%nspinup,   &
                                            GFS_Control%nlives_g,  &
                                            GFS_Control%ca_amplitude, &
                                            CA_grid)
              call remap_ca_to_fv3(CA_grid, GFS_Grid%xlon, GFS_Grid%xlat, fv3_ca)
              ca1_cpl = real(fv3_ca(:,:,1), kind_phys)
              ca2_cpl = 0.0
              ca3_cpl = 0.0
              if (GFS_Control%nca_g >= 2) ca2_cpl = real(fv3_ca(:,:,2), kind_phys)
              if (GFS_Control%nca_g >= 3) ca3_cpl = real(fv3_ca(:,:,3), kind_phys)

              if (GFS_Control%me == GFS_Control%master) then
                  print *, "CA: Outputting cellular automata"
                  call write_cellular_automata('ca_restart.nc', CA_grid, GFS_Control%kdt)
              endif
              call fv3ca_write(fv3_ca, GFS_Control%kdt, ierr)
              if (ierr /= 0) then
                call mpp_error(FATAL, 'Failed to write FV3 CA test output')
                return
              endif
          endif

!          call cellular_automata_global(GFS_Control%kdt,GFS_control%restart,GFS_Control%first_time_step,ca1_cpl,ca2_cpl,ca3_cpl,                &
!            Atm(mygrid)%domain_for_coupler, nblks,Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy,levs, &
!            GFS_Control%nca_g,GFS_Control%ncells_g,GFS_Control%nlives_g,GFS_Control%nfracseed,GFS_Control%nseed_g,                              &
!            GFS_Control%iseed_ca,GFS_control%tile_num,GFS_Control%ca_smooth,GFS_Control%nspinup,Atm_block%blksz(1),                             &
!            GFS_Control%nsmooth,GFS_Control%ca_amplitude,GFS_Control%master,GFS_Control%communicator)
          ! Copy contiguous data back
          do nb=1,nblks
             GFS_Coupling%ca1(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca1_cpl(nb,1:GFS_Control%blksz(nb))
             GFS_Coupling%ca2(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca2_cpl(nb,1:GFS_Control%blksz(nb))
             GFS_Coupling%ca3(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca3_cpl(nb,1:GFS_Control%blksz(nb))
          enddo
       endif

      endif !do_ca

    endif initalize_stochastic_physics

  contains
    
    ! Remaps gaussian grid to fv3 using stochy_la2ga
    subroutine remap_ca_to_fv3(gaussian_grid, target_lons, target_lats, fv3_fields)
      type(ca_grid_t), intent(in) :: gaussian_grid
      real(kind=kind_phys), intent(in) :: target_lons(:), target_lats(:)
      real(kind=kind_io8), intent(out) :: fv3_fields(:,:,:)
      real(kind=kind_io8), parameter :: rad2deg = 180.0_kind_io8 / acos(-1.0_kind_io8)
      integer :: block, field_index, blocksize, first_point, last_point

      fv3_fields = 0.0_kind_io8
      do field_index=1,gaussian_grid%nca
        do block=1,nblks
          blocksize = GFS_Control%blksz(block)
          first_point = GFS_Control%chunk_begin(block)
          last_point = GFS_Control%chunk_end(block)
          call stochy_la2ga(gaussian_grid%field(:,:,field_index), &
               gaussian_grid%nlon, gaussian_grid%nlat, gaussian_grid%lons, gaussian_grid%lats, &
               gaussian_grid%wlon, gaussian_grid%rnlat, fv3_fields(block,1:blocksize,field_index), &
               blocksize, real(target_lats(first_point:last_point), kind_io8) * rad2deg, &
               real(target_lons(first_point:last_point), kind_io8) * rad2deg)
        enddo
      enddo
    end subroutine remap_ca_to_fv3

    !! CA fv3 output for testing, can be removed!
    subroutine fv3ca_write(ca, time, ierr)
      use netcdf
      real(kind=kind_io8), intent(in) :: ca(:,:,:)
      integer, intent(in) :: time
      integer, intent(out) :: ierr
      character(len=128) :: filename
      character(len=32) :: var_name
      integer :: ncid, xdim_id, ydim_id, time_dim_id
      integer :: xt_var_id, yt_var_id, time_var_id
      integer :: nx_out, ny_out, nc_out, field_index, record, close_error
      integer :: block, point, grid_i, grid_j
      integer :: ca_ids(size(ca,3))
      real(kind=4) :: grid_xt(Atm_block%iec-Atm_block%isc+1)
      real(kind=4) :: grid_yt(Atm_block%jec-Atm_block%jsc+1)
      real(kind=4) :: workg(Atm_block%iec-Atm_block%isc+1, Atm_block%jec-Atm_block%jsc+1)
      real(kind=4), parameter :: undef = 9.99e+20
      logical :: file_exists

      nx_out = size(workg,1)
      ny_out = size(workg,2)
      nc_out = size(ca,3)
      grid_xt = [(real(point, kind=4), point=Atm_block%isc,Atm_block%iec)]
      grid_yt = [(real(point, kind=4), point=Atm_block%jsc,Atm_block%jec)]
      write(filename, '("RESTART/ca_out.pe", I6.6, ".nc")') GFS_Control%me
      inquire(file=trim(filename), exist=file_exists)

      if (.not. file_exists) then
        ierr = nf90_create(trim(filename), NF90_CLOBBER, ncid)
        if (ierr /= nf90_noerr) return
        ierr = nf90_def_dim(ncid, 'grid_xt', nx_out, xdim_id)
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_def_dim(ncid, 'grid_yt', ny_out, ydim_id)
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dim_id)
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_def_var(ncid, 'grid_xt', NF90_FLOAT, (/xdim_id/), xt_var_id)
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_put_att(ncid, xt_var_id, 'cartesian_axis', 'X')
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_def_var(ncid, 'grid_yt', NF90_FLOAT, (/ydim_id/), yt_var_id)
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_put_att(ncid, yt_var_id, 'cartesian_axis', 'Y')
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_def_var(ncid, 'time', NF90_INT, (/time_dim_id/), time_var_id)
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_put_att(ncid, time_var_id, 'units', 'model timesteps')
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_put_att(ncid, time_var_id, 'cartesian_axis', 'T')
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_put_att(ncid, NF90_GLOBAL, 'tile', GFS_Control%tile_num)
        if (ierr /= nf90_noerr) goto 900
        do field_index=1,nc_out
          write(var_name, '("ca", I0)') field_index
          ierr = nf90_def_var(ncid, trim(var_name), NF90_FLOAT, &
                             (/xdim_id,ydim_id,time_dim_id/), ca_ids(field_index))
          if (ierr /= nf90_noerr) goto 900
          ierr = nf90_put_att(ncid, ca_ids(field_index), 'long_name', 'random pattern' // trim(var_name))
          if (ierr /= nf90_noerr) goto 900
          ierr = nf90_put_att(ncid, ca_ids(field_index), 'units', 'None')
          if (ierr /= nf90_noerr) goto 900
          ierr = nf90_put_att(ncid, ca_ids(field_index), 'missing_value', undef)
          if (ierr /= nf90_noerr) goto 900
          ierr = nf90_put_att(ncid, ca_ids(field_index), '_FillValue', undef)
          if (ierr /= nf90_noerr) goto 900
          ierr = nf90_put_att(ncid, ca_ids(field_index), 'cell_methods', 'time: point')
          if (ierr /= nf90_noerr) goto 900
        enddo
        ierr = nf90_enddef(ncid)
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_put_var(ncid, xt_var_id, grid_xt)
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_put_var(ncid, yt_var_id, grid_yt)
        if (ierr /= nf90_noerr) goto 900
      else
        ierr = nf90_open(trim(filename), NF90_WRITE, ncid)
        if (ierr /= nf90_noerr) return
        ierr = nf90_inq_dimid(ncid, 'time', time_dim_id)
        if (ierr /= nf90_noerr) goto 900
        ierr = nf90_inq_varid(ncid, 'time', time_var_id)
        if (ierr /= nf90_noerr) goto 900
        do field_index=1,nc_out
          write(var_name, '("ca", I0)') field_index
          ierr = nf90_inq_varid(ncid, trim(var_name), ca_ids(field_index))
          if (ierr /= nf90_noerr) goto 900
        enddo
      endif

      ierr = nf90_inquire_dimension(ncid, time_dim_id, len=record)
      if (ierr /= nf90_noerr) goto 900
      record = record + 1
      ierr = nf90_put_var(ncid, time_var_id, time, start=(/record/))
      if (ierr /= nf90_noerr) goto 900
      do field_index=1,nc_out
        workg = undef
        do block=1,nblks
          do point=1,GFS_Control%blksz(block)
            grid_i = Atm_block%index(block)%ii(point) - Atm_block%isc + 1
            grid_j = Atm_block%index(block)%jj(point) - Atm_block%jsc + 1
            workg(grid_i,grid_j) = real(ca(block,point,field_index), kind=4)
          enddo
        enddo
        ierr = nf90_put_var(ncid, ca_ids(field_index), workg, start=(/1,1,record/))
        if (ierr /= nf90_noerr) goto 900
      enddo
900   close_error = nf90_close(ncid)
      if (ierr == nf90_noerr) ierr = close_error
    end subroutine fv3ca_write

    !> Transfers a 2D field to the stochastics module.
    !>
    !> This subroutine takes an input field and transfers it to the 
    !> stochastics module for further processing. The input data is 
    !> provided in `data_in` and the processed data is returned in 
    !> `data_out`.
    !>
    !> @param[in] blksz    The block size of the data.
    !> @param[in] data_in  The input 2D data field to be transferred.
    !> @param[out] data_out The output 2D data field after processing.
    !>
    !> @author Dominikus Heinzeller @date 8/6/21
    subroutine transfer_field_to_stochastics(blksz, data_in, data_out)

      integer, dimension(:), intent(in) :: blksz
      real(kind=kind_phys), dimension(:), intent(in) :: data_in
      real(kind=kind_phys), dimension(:,:), intent(out) :: data_out
      integer :: i, nb, ni

      nb = 1
      ni = 1
      do i=1,size(data_in)
        if (ni>blksz(nb)) then
          nb = nb+1
          ni = 1
        end if
        data_out(nb,ni) = data_in(i)
        ni =  ni+1
      end do

    end subroutine transfer_field_to_stochastics

    !> Transfers a 3D field to the stochastic physics module.
    !>
    !> This subroutine takes a 3D field as input and transfers it to the
    !> stochastic physics module for further processing. The input data
    !> is provided in `data_in` and the processed data is returned in `data_out`.
    !>
    !> @param[in] blksz The block size for the data transfer.
    !> @param[in] data_in The input 3D field data to be transferred.
    !> @param[out] data_out The output 3D field data after processing.
    !>
    !> @author Dominikus Heinzeller @date 8/6/21
    subroutine transfer_field_to_stochastics_3d(blksz, data_in, data_out)

      integer, dimension(:), intent(in) :: blksz
      real(kind=kind_phys), dimension(:,:), intent(in) :: data_in
      real(kind=kind_phys), dimension(:,:,:), intent(out) :: data_out
      integer :: j

      do j=1,size(data_in, dim=2)
         call transfer_field_to_stochastics(blksz, data_in(:,j), data_out(:,:,j))
      end do

    end subroutine transfer_field_to_stochastics_3d

    !> Transfers a 3D field from the stochastic physics module.
    !>
    !> @param[in] blksz The block size for the data transfer.
    !> @param[in] data_in The input 3D field data to be transferred.
    !> @param[out] data_out The output 3D field data after processing.
    !>  
    !> @author Dominikus Heinzeller @date 8/6/21
    subroutine transfer_field_from_stochastics(blksz, data_in, data_out)

      integer, dimension(:), intent(in) :: blksz
      real(kind=kind_phys), dimension(:,:), intent(in) :: data_in
      real(kind=kind_phys), dimension(:), intent(out) :: data_out
      integer :: i, nb, ni

      nb = 1
      ni = 1
      do i=1,size(data_out)
        if (ni>blksz(nb)) then
          nb = nb+1
          ni = 1
        end if
        data_out(i)= data_in(nb,ni)
        ni =  ni+1
      end do

    end subroutine transfer_field_from_stochastics

  end subroutine stochastic_physics_wrapper

  !> This subroutine finalizes the stochastic physics wrapper.
  !> It performs any necessary cleanup operations for the stochastic
  !> physics module.
  !>
  !> @param GFS_Control The control structure for the GFS model.
  !>
  !> @author Dominikus Heinzeller @date 8/6/21
  subroutine stochastic_physics_wrapper_end (GFS_Control)

  use GFS_typedefs,       only: GFS_control_type
  use stochastic_physics, only: finalize_stochastic_physics

  implicit none

  type(GFS_control_type),   intent(inout) :: GFS_Control

  if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type > 0) .OR. GFS_Control%do_spp) then
      if (allocated(xlat)) deallocate(xlat)
      if (allocated(xlon)) deallocate(xlon)
      if (GFS_Control%do_sppt) then
         if (allocated(sppt_wts)) deallocate(sppt_wts)
      end if
      if (GFS_Control%do_shum) then
         if (allocated(shum_wts)) deallocate(shum_wts)
      end if
      if (GFS_Control%do_skeb) then
         if (allocated(skebu_wts)) deallocate(skebu_wts)
         if (allocated(skebv_wts)) deallocate(skebv_wts)
      end if
      if (GFS_Control%do_spp) then
         if (allocated(spp_wts)) deallocate(spp_wts)
      end if
      if ( GFS_Control%lndp_type == 2 ) then
         lsoil = -999
         if (allocated(sfc_wts)) deallocate(sfc_wts)
      end if
      if (GFS_Control%lndp_type == 2) then
          if (allocated(smc))    deallocate(smc)
          if (allocated(slc))    deallocate(slc)
          if (allocated(stc))    deallocate(stc)
          if (allocated(stype))  deallocate(stype)
          if (allocated(vfrac))  deallocate(vfrac)
          if (allocated(snoalb)) deallocate(snoalb)
          if (allocated(alnsf))  deallocate(alnsf)
          if (allocated(alnwf))  deallocate(alnwf)
          if (allocated(semis))  deallocate(semis)
          if (allocated(zorll))  deallocate(zorll)
      endif
      call finalize_stochastic_physics()
   endif
   if(GFS_Control%do_ca)then
        if(GFS_Control%ca_sgs)then
           deallocate(sst         )
           deallocate(lmsk        )
           deallocate(lake        )
           deallocate(uwind       )
           deallocate(vwind       )
           deallocate(height      )
           deallocate(dx          )
           deallocate(condition   )
           deallocate(ca_deep_cpl )
           deallocate(ca_turb_cpl )
           deallocate(ca_shal_cpl )
        endif
        if(GFS_Control%ca_global)then
            deallocate(ca1_cpl )
            deallocate(ca2_cpl )
            deallocate(ca3_cpl )
        endif
   endif
  end subroutine stochastic_physics_wrapper_end

end module stochastic_physics_wrapper_mod
