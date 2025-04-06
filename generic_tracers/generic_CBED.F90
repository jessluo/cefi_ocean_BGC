module generic_CBED

  use g_tracer_utils, only : g_tracer_type, g_tracer_get_common
  use g_tracer_utils, only : g_tracer_set_values, g_tracer_get_values
  use g_tracer_utils, only : g_tracer_get_pointer
  use cobalt_types

implicit none ; private
! commment: checking if git is working
public generic_CBED_sediments_update_from_source

contains

  subroutine generic_CBED_sediments_update_from_source(tracer_list, cobalt, ilb, jlb, mask_coast, &
           grid_tmask, grid_dat, grid_kmt, isc,iec, jsc,jec, isd, jsd, nk, r_dt, dt, frunoff, rho_dzt, dzt, internal_heat)

    type(g_tracer_type),          pointer       :: tracer_list
    type(generic_COBALT_type),    intent(inout) :: cobalt
    type(phytoplankton), dimension(NUM_PHYTO), intent(inout) :: phyto
    integer,                      intent(in)    :: ilb, jlb
    real, dimension(ilb:,jlb:),   intent(in)    :: grid_dat
    real, dimension(:,:,:),       intent(in)    :: grid_tmask
    integer, dimension(:,:),      intent(in)    :: mask_coast, grid_kmt
    integer,                      intent(in)    :: isc,iec, jsc,jec, isd, jsd, nk
    real,                         intent(in)    :: r_dt, dt
    real, dimension(ilb:,jlb:),   intent(in)    :: frunoff
    real, dimension(ilb:,jlb:,:), intent(in)    :: rho_dzt, dzt
    real, dimension(ilb:,jlb:),   intent(in), optional :: internal_heat

    integer :: i, j, k
    real :: fpoc_btm, drho_dzt, log10_fpoc_btm
    integer, dimension(isc:iec,jsc:jec) :: k_bot
    real,    dimension(isc:iec,jsc:jec) :: rho_dzt_bot

	! CBED variables
	real, dimension(isc:iec,jsc:jec) :: fntot_in, fptot_in, ffetot_in, fsitot_in
	real, dimension(isc:iec,jsc:jec) :: no3, nh4, o2, dic, odu, alk
	real, dimension(isc:iec,jsc:jec) :: nomf, nomm, noms, nomtot   ! nitrogen in organic matter (nom) fast reacting, medium reacting, and slow reacting, and total
	real, dimension(isc:iec,jsc:jec) :: focf, focm, focs           ! flux of organic carbon fast, medium, slow (to calculate boundary) CHECK

  real, dimension(isc:iec,jsc:jec) :: o2resp_fast, o2resp_med, o2resp_slow        ! Rate of aerobic respiration 
  real, dimension(isc:iec,jsc:jec) :: no3resp_fast, no3resp_med, no3resp_slow     ! Rate of no3 respiration (denitrification)
  real, dimension(isc:iec,jsc:jec) :: anaerobicresp_fast, anaerobicresp_med, anaerobicresp_slow        ! Rate of anaerobic respiration 

  real, dimension(isc:iec,jsc:jec) :: nitrification, annamox, ODUox        ! secondary reactions 

	real, dimension(isc:iec,jsc:jec) :: gamma_fast, gamma_med, gamma_slow     ! organic matter decay rates, units of year-1
	
  real :: gamma_nitrif    =   1e6         ! rate constant for nitrification ( need to make sure the values are correct for the unit used )
  real :: gamma_anammox  =   1e6         ! rate constant for anammox
  real :: gamma_oduox    =   1e6         ! rate constant for ODU oxidation

  real :: k_o2  =   0.008            ! mol/m3        ! half-saturation constant for aerobic respiration (conc. unit)
  real :: k_no3 =   0.001            ! mol/m3        ! half-saturation constant for denitrification (conc. unit)

	real :: frac_omf = 0.7 ! fraction of organic matter in fast reacting pool
	real :: frac_omm = 0.2 ! fraction of organic matter in medium reacting pool
	real :: frac_oms = 0.1 ! fraction of organic matter in slow reacting pool
	real :: phi = 0.8      ! porosity, starting with a fixed value
	real :: Rho_solid = 2.5e3 ! solid density (kg/m3), use Rho_0 from cobalt for seawater density (1035 kg/m3)

    ! Calculate the bottom conditions and the fluxes to the bottom for diagnostics and benthic flux calculations.
    ! MOM4/5 used the bottom grid cell, but MOM6 often has a number of vanishingly thin layers overlying the bottom.
    ! Grid scale noise in these layers can occur, particularly for quantities with large bottom fluxes.  COBALT thus
    ! uses conditions over a specified bottom layer thickness (cobalt%bottom_thickness, default = 1m) for bottom calcs.

    do j = jsc, jec; do i = isc, iec  !{
       if (grid_kmt(i,j) .gt. 0) then !{

          ! Add the phytoplankton fluxes to the detritus fluxes to get total flux to benthos
          fntot_in(i,j) = cobalt%f_ndet_btf(i,j,1) + cobalt%f_ndi_btf(i,j,1) + &
            cobalt%f_nsm_btf(i,j,1) + cobalt%f_nmd_btf(i,j,1) + cobalt%f_nlg_btf(i,j,1)
          fptot_in(i,j) = cobalt%f_pdet_btf(i,j,1) + cobalt%f_pdi_btf(i,j,1) + &
            cobalt%f_psm_btf(i,j,1) + cobalt%f_pmd_btf(i,j,1) + cobalt%f_plg_btf(i,j,1)
          ffetot_in(i,j) = cobalt%f_fedet_btf(i,j,1) + cobalt%f_fedi_btf(i,j,1) + &
            cobalt%f_fesm_btf(i,j,1) + cobalt%f_femd_btf(i,j,1) + cobalt%f_felg_btf(i,j,1)
          fsitot_in(i,j) = cobalt%f_sidet_btf(i,j,1) + cobalt%f_silg_btf(i,j,1) + &
            cobalt%f_simd_btf(i,j,1)

		  ! Initialize tracers - TODO: start with initial conditions
		  nomf(i,j) = 0.0
		  nomm(i,j) = 0.0
		  noms(i,j) = 0.0
		  nomtot(i,j) = 0.0
		  
          ! Calculate the values of tracers influencing the sedimentary transformations
          ! and fluxes over a layer defined by "bottom_thickess".
          rho_dzt_bot(i,j) = 0.0
          o2(i,j) = 0.0
          no3(i,j) = 0.0
          nh4(i,j) = 0.0
          dic(i,j) = 0.0
          alk(i,j) = 0.0
		  
          k_bot(i,j) = 0
          ! Note that grid_kmt is always the total number of layers in MOM6
          do k = grid_kmt(i,j),1,-1   !{
            ! Check if the top of layer k is within the bottom thickness.  If so, include its properties in the bottom
            ! layer averages.  Overshoots will be subtracted off later.
            if (rho_dzt_bot(i,j).lt.(cobalt%Rho_0*cobalt%bottom_thickness)) then
              k_bot(i,j) = k
              rho_dzt_bot(i,j) = rho_dzt_bot(i,j) + rho_dzt(i,j,k)
              o2(i,j) = o2(i,j) + cobalt%f_o2(i,j,k)*rho_dzt(i,j,k) 
              no3(i,j) = no3(i,j) + cobalt%f_no3(i,j,k)*rho_dzt(i,j,k) 
              nh4(i,j) = nh4(i,j) + cobalt%f_nh4(i,j,k)*rho_dzt(i,j,k) 
              dic(i,j) = dic(i,j) + cobalt%f_dic(i,j,k)*rho_dzt(i,j,k) 
              alk(i,j) = alk(i,j) + cobalt%f_alk(i,j,k)*rho_dzt(i,j,k) 
            endif
          enddo
          ! Subtract off overshoot
          drho_dzt = rho_dzt_bot(i,j) - cobalt%Rho_0*cobalt%bottom_thickness
          o2(i,j)=o2(i,j)-cobalt%f_o2(i,j,k_bot(i,j))*drho_dzt
          no3(i,j)=no3(i,j)-cobalt%f_no3(i,j,k_bot(i,j))*drho_dzt
          nh4(i,j)=nh4(i,j)-cobalt%f_nh4(i,j,k_bot(i,j))*drho_dzt
          dic(i,j)=dic(i,j)-cobalt%f_dic(i,j,k_bot(i,j))*drho_dzt
          alk(i,j)=alk(i,j)-cobalt%f_alk(i,j,k_bot(i,j))*drho_dzt
          ! convert back to moles kg-1
		  o2(i,j)=o2(i,j)/(cobalt%bottom_thickness*cobalt%Rho_0)
          no3(i,j)=no3(i,j)/(cobalt%bottom_thickness*cobalt%Rho_0)
          nh4(i,j)=nh4(i,j)/(cobalt%bottom_thickness*cobalt%Rho_0)
          dic(i,j)=dic(i,j)/(cobalt%bottom_thickness*cobalt%Rho_0)
          alk(i,j)=alk(i,j)/(cobalt%bottom_thickness*cobalt%Rho_0)
		 

          ! Calculate the processing of organic matter in the sediment.  The fate of organic matter is partitioned
          ! between burial (i.e., removal from the system), aerobic remineralization, remineralization via 
          ! denitrification, and remineralization via sulfate reduction.  Note that the latter pathway is effectively
          ! a "catch all" for any other anaerobic pathway and the sulfate cycle is not explicitly modeled.
          k = grid_kmt(i,j)
          if (fntot_in(i,j) .gt. 0.0) then !{
			 
			 focf(i,j) = fntot_in(i,j)*cobalt%c_2_n * frac_omf * Rho_solid ! units of mol/m3 organic matter
			 focm(i,j) = fntot_in(i,j)*cobalt%c_2_n * frac_omm * Rho_solid
			 focs(i,j) = fntot_in(i,j)*cobalt%c_2_n * frac_oms * Rho_solid
			 
       ! calculate converstion from 'foc' to 'oc' conc in sediment  NEED TO DO 
			 ! Should the reactions be written in 'oc' or use 'on' with c_2_n ?

			 ! aerobic respiration
			 o2resp_fast(i,j) = gamma_fast * nomf(i,j) * o2(i,j)/(k_o2 + o2(i,j))
			 o2resp_med(i,j) = gamma_med * nomm(i,j) * o2(i,j)/(k_o2 + o2(i,j))
			 o2resp_slow(i,j) = gamma_slow * noms(i,j) * o2(i,j)/(k_o2 + o2(i,j))
			 
       ! denitrification 
       no3resp_fast(i,j)   = gamma_fast * nomf(i,j)  *  no3(i,j)/(k_no3 + no3(i,j)) * k_o2/(k_o2 + o2(i,j))
			 no3resp_med(i,j)    = gamma_med  * nomm(i,j)  *  no3(i,j)/(k_no3 + no3(i,j)) * k_o2/(k_o2 + o2(i,j))
			 no3resp_slow(i,j)   = gamma_slow * noms(i,j)  *  no3(i,j)/(k_no3 + no3(i,j)) * k_o2/(k_o2 + o2(i,j))

       ! anaerobic respiration  
       anaerobicresp_fast(i,j)  = gamma_fast * nomf(i,j)  * k_no3/(k_no3 + no3(i,j)) * k_o2/(k_o2 + o2(i,j))
			 anaerobicresp_med(i,j)   = gamma_med  * nomm(i,j)  * k_no3/(k_no3 + no3(i,j)) * k_o2/(k_o2 + o2(i,j))
			 anaerobicresp_slow(i,j)  = gamma_slow * noms(i,j)  * k_no3/(k_no3 + no3(i,j)) * k_o2/(k_o2 + o2(i,j))

       !! secondary reactions

       ! nitrification 
       nitrification(i,j) = gamma_nitrif * nh4(i,j) * o2(i,j)

       ! annamox 
       annamox(i,j) = gamma_anammox * nh4(i,j) * no3(i,j)

       ! ODU oxidation 
       ODUox(i,j) = gamma_oduox * odu(i,j) * o2(i,j)


          endif !}

          !
          ! Bottom flux boundaries passed to the vertical mixing routine
          ! (negative values are fluxes into the ocean)
          !
          cobalt%b_dic(i,j) = 0.0
          cobalt%b_fed(i,j) = 0.0
          cobalt%b_nh4(i,j) = 0.0
          cobalt%b_no3(i,j) = 0.0
          cobalt%b_o2(i,j)  = 0.0
          cobalt%b_alk(i,j) = 0.0
          cobalt%b_po4(i,j) = 0.0
          cobalt%b_sio4(i,j)= 0.0

       endif !}
    enddo; enddo  !} i, j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%f_cased(i,j,k) = 0.0
    enddo; enddo ; enddo  !} i,j,k

    call g_tracer_set_values(tracer_list,'alk',  'btf', cobalt%b_alk ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic',  'btf', cobalt%b_dic ,isd,jsd)
    call g_tracer_set_values(tracer_list,'fed',  'btf', cobalt%b_fed ,isd,jsd)
    call g_tracer_set_values(tracer_list,'nh4',  'btf', cobalt%b_nh4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'no3',  'btf', cobalt%b_no3 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'o2',   'btf', cobalt%b_o2  ,isd,jsd)
    call g_tracer_set_values(tracer_list,'po4',  'btf', cobalt%b_po4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'sio4', 'btf', cobalt%b_sio4,isd,jsd)

  end subroutine generic_CBED_sediments_update_from_source

end module generic_CBED
